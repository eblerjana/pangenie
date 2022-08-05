import sys
import argparse
import gzip
from variantclassifier import VariantType, determine_variant_from_line
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np 

def check_mendelian_consistency(child_gt, parent1_gt, parent2_gt):
	child = child_gt.get_alleles()
	parent1 = parent1_gt.get_alleles()
	parent2 = parent2_gt.get_alleles()
	if child[0] in parent1 and child[1] in parent2:
		return True
	if child[0] in parent2 and child[1] in parent1:
		return True
	return False

class Genotype:
	def __init__(self, alleles, is_phased):
		self._alleles = alleles
		self._phased = is_phased

	def __str__(self):
		if self._phased:
			return '|'.join([str(i) for i in self._alleles])
		else:
			return '/'.join([str(i) for i in self._alleles])

	def get_alleles(self):
		return self._alleles

	def get_ploidy(self):
		return len(self._alleles)

	def is_phased(self):
		return self._phased

	def __eq__(self, other):
		return sorted(self._alleles) == sorted(other._alleles)

	def is_hom_ref(self):
		return all([int(a) == 0 for a in self._alleles])

	def is_none(self):
		return self._alleles == []

def genotype_from_string(gt_string):
	is_phased = False
	alleles = []
	if '.' in gt_string:
		# untyped
		return Genotype(alleles, is_phased)

	if '|' in gt_string:
		is_phased = True
		alleles = [int(allele) for allele in gt_string.split('|')]
	elif '/' in gt_string:
		is_phased = False
		alleles = [int(allele) for allele in gt_string.split('/')]
	else:
		assert False
	return Genotype(alleles, is_phased)

def parse_trios(ped_file, samples):
	samples_to_include = set()
	with open(samples, 'r') as listing:
		for line in listing:
			samples_to_include.add(line.strip())
	trios = {}
	for line in open(ped_file, 'r'):
		if line.startswith('#'):
			continue
		fields = line.split()
		# map trio_name -> [child, parent1, parent2, nr consistent, nr_inconsistent, nr_untyped, nr total]
		if any([s not in samples_to_include for s in fields[1:4]]):
			continue
		trios[(fields[1], fields[2], fields[3])] = [fields[1], fields[2], fields[3], 0, 0, 0, 0]
	return trios


def remove_samples_from_header(header_line, samples):
	fields = header_line.split()
	result = fields[:9]
	for s in fields[9:]:
		if s not in samples:
			result.append(s)
	return '\t'.join(result)

def check_all_absent(genotype_strings):
	# TODO this assumes that the GT is the first field in sample column
	genotypes = [genotype_from_string(g.split(':')[0]) for g in genotype_strings]
	skip_line=True
	for genotype in genotypes:
		if not genotype.is_hom_ref():
			skip_line = False
			break
	return skip_line

def run_statistics(vcf_file, ped_file, samples, table_tsv, column_prefix):
	# map child -> [parents]
	trios = parse_trios(ped_file, samples)
	# histograms[vartype][i] = number of variants consistent in i trios
	histograms = { vartype.name : [0]*(len(trios)+1) for vartype in VariantType }
	histograms['multi_var'] = [0] * (len(trios)+1)
	untyped_vars = 0
	multitype_vars = 0

	with open(table_tsv, 'w') as out_tsv:
		header = '\t'.join(['variant_id', column_prefix + '_mendelian_consistent_trios', column_prefix + '_alternative_transmitted', column_prefix + '_considered_trios', column_prefix + '_all_0/0', column_prefix + '_all_0/1', column_prefix + '_all_1/1', column_prefix + '_untyped_alleles_present'])
		out_tsv.write(header + '\n')
		header = None
		sample_to_index = {}
		for record in gzip.open(vcf_file, 'rt'):
			if record.startswith('##'):
				continue
			if record.startswith('#'):
				header = record.strip().split()
				for i,f in enumerate(header):
					sample_to_index[f] = i
				continue
			assert header is not None
			fields = record.split()
			genotype_index = fields[8].split(':').index('GT')
			consistent_trios = 0
			total_consistent_trios = 0
			alt_transmitted = 0
			all_het = 0
			all_abs = 0
			all_present = 0
			assert len(fields[4].split(',')) == 1
			info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i} 
			assert 'ID' in info_fields
			assert len(info_fields['ID'].split(',')) == 1
			variant_id = info_fields['ID']
			total_trios = 0
			untyped_trios = 0
			if 'X' in fields[0] or 'Y' in fields[0]:
				continue
			vartype = determine_variant_from_line(record)
			assert isinstance(vartype, VariantType), 'Unexpected return value: {}'.format(type(vartype))
			untyped = False
			for name in trios:
				field_child = fields[sample_to_index[trios[name][0]]].split(':')
				field_parent1 = fields[sample_to_index[trios[name][1]]].split(':')
				field_parent2 = fields[sample_to_index[trios[name][2]]].split(':')
				gt_child = genotype_from_string(field_child[genotype_index])
				gt_parent1 = genotype_from_string(field_parent1[genotype_index])
				gt_parent2 = genotype_from_string(field_parent2[genotype_index])
				if any([g.is_none() for g in [gt_child, gt_parent1, gt_parent2]]):
					untyped_trios += 1
					if not untyped:
						untyped = True
						untyped_vars += 1
				elif all([ g == gt_child for g in [gt_parent1, gt_parent2]]):
					# all genotypes same, automatically consistent
					if gt_child == Genotype([0,0], False):
						all_abs += 1
					elif gt_child == Genotype([0,1], False):
						all_het += 1
					elif gt_child == Genotype([1,1], False):
						all_present += 1
					else:
						assert(False)
					total_consistent_trios += 1
				else:
					total_trios += 1
					consistent = check_mendelian_consistency(gt_child, gt_parent1, gt_parent2)
					if consistent:
						consistent_trios += 1
						total_consistent_trios += 1
						# check how often alt allele was transmitted
						alt_transmitted += sum(a!=0 for a in gt_child.get_alleles())
			assert total_trios + all_abs + all_het + all_present + untyped_trios == len(trios)
			if not untyped:
				histograms[vartype.name][total_consistent_trios] += 1
			# write to output
			out_tsv.write('\t'.join([variant_id, str(consistent_trios), str(alt_transmitted), str(total_trios), str(all_abs), str(all_het), str(all_present), str(untyped_trios)]) + '\n')
	# output the results
	print('variants_untyped\t' + str(untyped_vars))
#	print('variants_multitype\t' + str(multitype_vars))
	print('\t'.join(['nr_trios'] + [str(i) for i in range(0,len(trios)+1)]))
	for vartype in VariantType:
		print('\t'.join([vartype.name] + [str(i) for i in histograms[vartype.name]]))
		

def run_mendelian(vcf, ped, samples, output, remove_children=False):
	chromosomes = ['chr' + str(i) for i in range(1,23)] + ['chrX']
	# map name -> [trio]
	trios = parse_trios(ped, samples)
	# map sample -> column in vcf
	sample_to_id = None
	child_ids = None
	for line in open(vcf, 'r'):
		if line.startswith("##"):
			print(line[:-1])
			continue
		if line.startswith("#"):
			sample_to_id = {}
			fields = line.split()
			assert len(fields) > 9
			for i, sample in enumerate(fields[9:]):
				sample_to_id[sample] = i + 9
			children = [i[0] for i in trios.values()]
			child_ids = [sample_to_id[i] for i in children]
			# write header line
			samples_to_skip = children if remove_children else []
			header_line = remove_samples_from_header(line, samples_to_skip)
			print(header_line)
			continue
		assert sample_to_id is not None
		fields = line.split()
		if not fields[0] in chromosomes:
			continue
		sex_chrom = False
		if 'X' in fields[0] or 'Y' in fields[0]:
			sex_chrom = True
#			continue
		# check mendelian consistency for all trios
		all_consistent = True
		for name in trios:
			if sex_chrom:
				continue
			child = trios[name][0]
			parent1 = trios[name][1]
			parent2 = trios[name][2]
			# determine genotypes
			# TODO this assumes that the GT is the first field in sample column
			child_gt = genotype_from_string(fields[sample_to_id[child]].split(':')[0])
			p1_gt = genotype_from_string(fields[sample_to_id[parent1]].split(':')[0])
			p2_gt = genotype_from_string(fields[sample_to_id[parent2]].split(':')[0])
			untyped = any(g.is_none() for g in [child_gt,p1_gt,p2_gt])
			# check mendelian consistency [parent1, parent2, nr consistent, nr_inconsistent, nr_untyped, nr total]
			if untyped:
				consistent = False
				trios[name][5] += 1 
			else:
				consistent = check_mendelian_consistency(child_gt, p1_gt, p2_gt)
			if not consistent:
				if not untyped: 
					trios[name][4] += 1
					sys.stderr.write(' '.join(['mendelian conflict in ', '-'.join(name), ' at position', fields[0], fields[1]+'\n']))
				all_consistent = False
			trios[name][6] += 1
			if consistent:
				trios[name][3] += 1
		if all_consistent or sex_chrom:
			only_parents = [f for i,f in enumerate(fields) if not i in child_ids]
			# variants only in children where removed as they cause mendelian conflicts
			if not check_all_absent(fields[9:]):
				assert (not check_all_absent(only_parents[9:])) or sex_chrom
			if remove_children:
				print('\t'.join(only_parents))
			else:
				print('\t'.join(fields))

	
	# write output
	with open(output, 'w') as outfile:
		outfile.write('\t'.join(['#sample', 'consistent_genotypes', 'inconsistent_genotypes', 'untyped', 'total_genotypes', 'mendelian_consistency\n']))
		for trio in trios:
			consistent = float(trios[trio][3])
			inconsistent = float(trios[trio][4])
			untyped = float(trios[trio][5])
			total = float(trios[trio][6])
			assert consistent+inconsistent+untyped == total
			mendelian_consistency = consistent / total 
			outfile.write('\t'.join(['-'.join(trio), str(consistent), str(inconsistent), str(untyped), str(total), str(mendelian_consistency)]) + '\n')


def run_plot_filter(outprefix, tsv):
	consistencies = []
	names = []
	for line in open(tsv, 'r'):
		if line.startswith('#'):
			continue
		fields = line.split()
		names.append(fields[0])
		consistencies.append(float(fields[-1])*100.0)
	x_values = [i for i in range(len(consistencies))]
	plt.figure(figsize=(15, 4.8))
	plt.bar(x_values, consistencies, width=0.9)
	plt.xticks(x_values, names, rotation=90)
	plt.ylabel('mendelian consistency [%]')
	plt.tight_layout()
	plt.savefig(outprefix + '_filter.pdf')
	plt.close()

def run_plot_statistics(outprefix, tsv):
	type_to_numbers = {}
	nr_trios = -1
	for line in open(tsv, 'r'):
		if line.startswith('variants_untyped') or line.startswith('nr_trios') or line.startswith('variants_multitype'):
			continue
		fields = line.split()
		nr_trios = len(fields)-1
		type_to_numbers[fields[0]] = [int(i) for i in fields[1:]]
	all_types = [0] * nr_trios
	x_values = [i for i in range(nr_trios)]
	for var, numbers in type_to_numbers.items():
		assert len(numbers) == nr_trios
		plt.bar(x_values, numbers)
		plt.title(var)
		plt.ylabel('count')
		plt.yscale('log')
		plt.xlabel('nr of mendelian consistent trios')
#		plt.xticks(x_values)
		plt.savefig(outprefix + '_statistics_' + var + '.pdf')
		plt.close()
		for i,number in enumerate(numbers):
			all_types[i] += number
	#print(type_to_numbers)
	plt.bar(x_values, all_types)
	#print(all_types, all_types[-1]/sum(all_types))
	plt.title('all variants')
	plt.ylabel('count')
	plt.yscale('log')
	plt.xlabel('nr of mendelian consistent trios')
	plt.xticks(x_values)
	plt.savefig(outprefix + '_statistics_all.pdf')
	plt.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='mendelian-consistency.py', description=__doc__)
	subparsers = parser.add_subparsers(dest="subparser_name")

	parser_mendelian = subparsers.add_parser('filter', help='check mendelian consistency for specified trios and remove inconsistent variants.')
	parser_mendelian.add_argument('-vcf', metavar='VCF', required=True, help='Multisample VCF-file with genotypes.')
	parser_mendelian.add_argument('-ped', metavar='PED', required=True, help='Trio relationships to be considered.')
	parser_mendelian.add_argument('-samples', metavar='SAMPLES', required=True, help='Samples to include')
	parser_mendelian.add_argument('-o', metavar='STATISTICS', required=True, help='File to write statistics to.')
	parser_mendelian.add_argument('--remove-children', default=False, action='store_true', help='Skip children in final VCF.')

	parser_statistics = subparsers.add_parser('statistics', help='print statistics.')
	parser_statistics.add_argument('-vcf', metavar='VCF', required=True, help='Multisample VCF-file with genotypes.')
	parser_statistics.add_argument('-ped', metavar='PED', required=True, help='Trio relationships to be considered.')
	parser_statistics.add_argument('-samples', metavar='SAMPLES', required=True, help='Samples to include')
	parser_statistics.add_argument('-table', metavar='TABLE', required=True, help='Write statistics per variant ID.')
	parser_statistics.add_argument('-column-prefix', metavar='COLUMNPREFIX', required=True, help='Column prefix in output.')

	parser_plot = subparsers.add_parser('plot', help='create plots from outputs of commands statistics and filter.')
	parser_plot.add_argument('-statistics', metavar='TSV', help='tsv produced by command statistics.')
	parser_plot.add_argument('-filter', metavar='TSV', help='tsv produced by command filter.')
	parser_plot.add_argument('outprefix', metavar='OUTPREFIX', help='prefix of the output files.')

	args = parser.parse_args()

	if args.subparser_name == 'filter':
		run_mendelian(args.vcf, args.ped, args.samples, args.o, args.remove_children)
	if args.subparser_name == 'statistics':
		run_statistics(args.vcf, args.ped, args.samples, args.table, args.column_prefix)
	if args.subparser_name == 'plot':
		if args.filter:
			run_plot_filter(args.outprefix, args.filter)
		if args.statistics:
			run_plot_statistics(args.outprefix, args.statistics)
