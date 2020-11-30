#!/usr/bin/python

import sys
from collections import defaultdict
from collections import namedtuple
import argparse
import math
from decimal import Decimal 
import vcf
from enum import Enum

Position = namedtuple("Position", "chrom position")

class VariantType(Enum):
	snp = 0
	small_insertion = 1
	small_deletion = 2
	small_complex = 3
	midsize_insertion = 4
	midsize_deletion = 5
	midsize_complex = 6
	large_insertion = 7
	large_deletion = 8
	large_complex = 9


class Variant:
	def __init__(self, binary_genotype, quality=None):
		"""
		Represents a genotyped variant.
		"""
		self._binary_genotype=binary_genotype
		self._quality=quality

	def get_binary_genotype(self):
		return self._binary_genotype

	def get_quality(self):
		return self._quality


def determine_type_from_ids(ids):
	"""
	Determine variant type from variant ID.
	"""
	ids = set(ids)
	results = []
	for var in ids:
		results.append(id_to_vartype(var))	
	return results


def id_to_vartype(var):
	var_type = var.split('-')[2]
	var_len = 1 if var_type == "SNV" else int(var.split('-')[-1])
	if var_type == "SNV":
		return VariantType.snp
	if var_len < 20:
		if var_type == "INS":
			return VariantType.small_insertion
		elif var_type == "DEL":
			return VariantType.small_deletion
		else:
			return VariantType.small_complex
	if var_len >= 20 and var_len < 50:
		if var_type == "INS":
			return VariantType.midsize_insertion
		elif var_type == "DEL":
			return VariantType.midsize_deletion
		else:
			return VariantType.midsize_complex
	if var_len >= 50:
		if var_type == "INS":
			return VariantType.large_insertion
		elif var_type == "DEL":
			return VariantType.large_deletion
		else:
			return VariantType.large_complex


def determine_genotypes_from_ids(ids, gt):
	allele_to_variants = {0:[]}
	variants_to_genotype = {}
	for i,var_id in enumerate(ids):
		allele_to_variants[i+1] = []
		for single_id in var_id.split(':'):
			allele_to_variants[i+1].append( (single_id, id_to_vartype(single_id)) )
			variants_to_genotype[single_id,id_to_vartype(single_id)] = 0

	genotype_string = gt.replace('/', '|')
	genotype_list = genotype_string.split('|')
	assert 1 <= len(genotype_list) <= 2
	if len(genotype_list) == 1:
		# haploid genotype
		genotype_list.append(genotype_list[0])
	for allele in genotype_list:
		if allele == '.':
			continue
		for (id,type) in allele_to_variants[int(allele)]:
			variants_to_genotype[(id,type)] += 1
			assert variants_to_genotype[(id,type)] < 3
	return variants_to_genotype


def extract_call(line, vcf_samples, samples, read_gl=False, read_qual=False):
	"""
	Extract genotype information for all samples from a VCF Record.
	"""
	# store var_ID -> sample -> genotype
	result = defaultdict(list)
	# vcf fields
	fields = line.split()
	# info fields present in vcf
	info_fields = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';')}
	# get all variant IDs/types	
	assert 'ID' in info_fields
	allele_ids = set([])
	for i in info_fields['ID'].split(','):
		for j in i.split(':'):
			allele_ids.add(j)
	# parse sample columns and store genotype info for each sample
	format_fields = fields[8].split(':')
	sample_columns = defaultdict(dict)
	for i,vcf_sample in enumerate(vcf_samples):
		genotype_info = fields[9+i].split(':')
		for f,v in zip(format_fields,genotype_info):
			sample_columns[vcf_sample][f] = v	
	# sample specific information
	for sample in samples:
		# genotype quality
		likelihood = 0
		genotype = sample_columns[sample]['GT']
		if (genotype[0] != '.') and (genotype[-1] != '.'):
			# compute likelihood of genotype
			if read_gl:
				# check if GQ field is present
				if 'GQ' in format_fields:
					likelihood = int(sample_columns[sample]['GQ'])
			if read_qual:
				if fields[5] != '.':
					likelihood = int(fields[5])
			# extract genotype information per variant ID
			all_genotypes = determine_genotypes_from_ids(info_fields['ID'].split(','), sample_columns[sample]['GT'])
			assert len(allele_ids) == len(all_genotypes)
			for (var_id,var_type) in all_genotypes:
				result[var_id].append(Variant(all_genotypes[(var_id, var_type)], likelihood))
		else:
			for i in allele_ids:
				result[i].append(Variant(3, likelihood))
	return result


class GenotypingStatistics:
	def __init__(self, var_id, quality):
		"""
		Collect some statistics about the genotypes.
		"""
		self.var_id = var_id
		# thresholds used to filter variant set
		self.quality = quality
		# numbers of variants in baseline and callset
		self.total_baseline = 0
		self.total_intersection = 0
		# fractions of correct, wrong and untyped variants
		self.correct_all = 0
		self.wrong_all = 0
		self.not_typed_all = 0
		# confusion matrix
		self.confusion_matrix = [[0 for x in range(4)] for y in range(4)]
		# absent in ground truth
		self.absent_truth = 0
		self.not_in_callset_all = 0


	def print_variant_stats_to_file(self, tsv_file):
		typed_all = max(self.correct_all + self.wrong_all, 1)
		assert self.total_baseline == self.correct_all + self.wrong_all + self.not_typed_all + self.not_in_callset_all
		# linear representation of the confusion matrix
		matrix = []
		for true in range(3):
			for typed in range(4):
				matrix.append(str(self.confusion_matrix[true][typed]))
		tsv_file.write('\t'.join([	self.var_id,
						str((self.correct_all/float(typed_all))*100.0), # correct
						str((self.wrong_all/float(typed_all))*100.0), # wrong
						str(((self.not_typed_all)/float(self.total_baseline if self.total_baseline is not 0 else 1))*100.0), # not typed
						str(self.correct_all),
						str(self.wrong_all),
						str(self.not_typed_all),
						str(self.absent_truth)	] + matrix) + '\n'
						)

	def print_matrix_to_file(self, txt_file):
		txt_file.write('\n###############################################################################################################################################\n')
		txt_file.write('#       Matrix for variant: '+ self.var_id + ' quality: '  + str(self.quality) + '   #\n')
		txt_file.write('###############################################################################################################################################\n\n')
		matrix_string = "\t0\t1\t2\t./.\n"
		for true in range(4):
			line = str(true) + '\t'
			for geno in range(4):
				line += str(self.confusion_matrix[true][geno]) + '\t'
			matrix_string += line + '\n'
		txt_file.write(matrix_string)


class GenotypeConcordanceComputer:
	def __init__(self, baseline_vcf, callset_vcf, samples, use_qual):
		"""
		Read variants from baseline and callset VCFs and store them
		for evaluation.
		"""
		self._baseline_variants = {}
		self._callset_variants = defaultdict(list)
		self._total_baseline = 0
		self._total_callset = 0
		self._samples = samples

		# read baseline variants
		baseline_samples = None
		with open(baseline_vcf) as baseline_vars:
			for line in baseline_vars:
				if line.startswith('##'):
					# description line, skip
					continue
				if line.startswith('#'):
					# header line, get list of samples 
					baseline_samples = line.split()[9:]
					continue
				assert baseline_samples is not None
				# determine genotypes per variant ID
				ids = extract_call(line, baseline_samples, self._samples)
				for var_id, genotypes in ids.items():
					self._baseline_variants[var_id] = genotypes
				self._total_baseline += len(ids)
		# read callset variants
		callset_samples = None
		with open(callset_vcf) as callset_vars:
			for line in callset_vars:
				if line.startswith('##'):
					# description line, skip
					continue
				if line.startswith('#'):
					# header line, get list of samples
					callset_samples = line.split()[9:]
					continue
				assert callset_samples is not None
				if use_qual:
					ids = extract_call(line, callset_samples, self._samples, read_qual = True)
				else:
					ids = extract_call(line, callset_samples, self._samples, read_gl = True)
				for var_id, genotypes in ids.items():
					self._callset_variants[var_id] = genotypes
				self._total_callset += len(ids)
		sys.stderr.write('Read ' + str(self._total_baseline) + ' variants from the baseline VCF.\n')
		sys.stderr.write('Read ' + str(self._total_callset) + ' variants from the callset VCF.\n')

		print(sys.getsizeof(self._callset_variants), sys.getsizeof(self._baseline_variants))

	def get_samples(self):
		return self._samples
	
	def print_variant_statistics(self, quality, tsv_file, txt_file):
		"""
		Compute genotype concordance statistics for each single variant across all samples.
		"""	
		for var_id, sample_to_genotype in self._baseline_variants.items():
			var_type = id_to_vartype(var_id)
			statistics = GenotypingStatistics(var_id, quality)
			callset_predictions = self._callset_variants[var_id]
#			if len(callset_predictions) > 1:
#				print('Warning: multiple predictions for variant ID ' + var_id + '. Using first.')
			for sample_id,sample in enumerate(self._samples):
				baseline_gt = sample_to_genotype[sample_id]
				if baseline_gt.get_binary_genotype() == 3:
					# ground truth unknown, skip variant
					statistics.absent_truth += 1
					continue
				statistics.total_baseline += 1
				# check genotype predictions of each sample
				if len(callset_predictions) > 0:
					callset_gt = callset_predictions[sample_id]
					statistics.total_intersection += 1
					if callset_gt.get_quality() < quality:
						print('untyped', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
						statistics.not_typed_all += 1
						continue
					statistics.confusion_matrix[baseline_gt.get_binary_genotype()][callset_gt.get_binary_genotype()] += 1
					# untyped genotype
					if callset_gt.get_binary_genotype() == 3:
						statistics.not_typed_all += 1
						print('untyped', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
						continue
					
					# check if genotype predictions are identical
					if callset_gt.get_binary_genotype() == baseline_gt.get_binary_genotype():
						statistics.correct_all += 1
						print('correct', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
					else:
						print('wrong', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
						statistics.wrong_all += 1				
				else:
					statistics.not_in_callset_all += 1
					statistics.confusion_matrix[baseline_gt.get_binary_genotype()][3] += 1
			statistics.print_variant_stats_to_file(tsv_file)
			statistics.print_matrix_to_file(txt_file)


def compute_stats_per_variant(baseline, callset, samples, use_qual, file_prefix, column_prefix, quality_thresholds):
	"""
	Compute variant specific statistics for given variant types.
	"""
	genotype_concordance = GenotypeConcordanceComputer(args.baseline, args.callset, samples, args.use_qual)
	print(sys.getsizeof(genotype_concordance))
	# sample specific evaluation
	tsv_file = open(file_prefix + '_variant-stats.tsv', 'w')
	txt_file = open(file_prefix + '_variant-stats.txt', 'w')
	# write headers to tsv files
	header = '\t'.join([	'variant_id',
				column_prefix + '_correct [%]',
				column_prefix + '_wrong [%]',
				column_prefix + '_not_typed [%]',
				column_prefix + '_correct',
				column_prefix + '_wrong',
				column_prefix + '_not_typed',
				column_prefix + '_absent_in_truth',
				column_prefix + '_0/0_typed_0/0',
				column_prefix + '_0/0_typed_0/1',
				column_prefix + '_0/0_typed_1/1',
				column_prefix + '_0/0_typed_./.',
				column_prefix + '_0/1_typed_0/0',
				column_prefix + '_0/1_typed_0/1',
				column_prefix + '_0/1_typed_1/1',
				column_prefix + '_0/1_typed_./.',
				column_prefix + '_1/1_typed_0/0',
				column_prefix + '_1/1_typed_0/1',
				column_prefix + '_1/1_typed_1/1',
				column_prefix + '_1/1_typed_./.'
							]) + '\n'
	tsv_file.write(header)
	genotype_concordance.print_variant_statistics(0, tsv_file, txt_file)
	tsv_file.close()
	txt_file.close()

if __name__ == "__main__":

	# baseline: multisample VCF containing variants and panel of true genotypes
	# callset: multisample VCF containing genotype predictions per sample
	parser = argparse.ArgumentParser(prog='genotype-concordance-variant.py', description=__doc__)
	parser.add_argument('baseline', metavar='BASELINE', help='multisample baseline VCF (ground truth genotypes).')
	parser.add_argument('callset', metavar='CALLSET', help='multisample callset VCF (genotyped variants).')
	parser.add_argument('file_prefix', metavar='OUTFILE', help='prefix of the output file name.')
	parser.add_argument('samples', metavar='SAMPLES', help='comma separated list of samples to evaluate.')
	parser.add_argument('column_prefix', metavar='PREFIX', help='prefix of the output column names.')
	parser.add_argument('--qualities', default='0', metavar='GQ-THRESHOLDS', help='comma separated list of GQ-thresholds to consider (default: consider all variants regardless of quality).')
	parser.add_argument('--use-qual', default=False, action='store_true', help='use qualities in QUAL field instead of GQ fields (default: use GQ field).')
	args = parser.parse_args()

	quality_thresholds = [int(i) for i in args.qualities.split(',')]
	samples = args.samples.split(',')	
	compute_stats_per_variant(args.baseline, args.callset, samples, args.use_qual, args.file_prefix, args.column_prefix, quality_thresholds)
