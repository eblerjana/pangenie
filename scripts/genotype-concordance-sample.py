#!/usr/bin/python

import sys
from collections import defaultdict
from collections import namedtuple
import argparse
import math
from decimal import Decimal 
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
	def __init__(self, binary_genotype, quality=None, allele_frequency=None, unique_kmers=None, missing_alleles=None):
		"""
		Represents a genotyped variant.
		"""
		self._binary_genotype=binary_genotype
		self._quality=quality
		self._allele_frequency=allele_frequency
		self._unique_kmers=unique_kmers
		self._missing_alleles=missing_alleles

	def get_binary_genotype(self):
		return self._binary_genotype

	def get_quality(self):
		return self._quality

	def get_allele_frequency(self):
		return self._allele_frequency

	def get_unique_kmers(self):
		return self._unique_kmers

	def get_missing_alleles(self):
		return self._missing_alleles

def determine_type_from_ids(ids):
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

	# allele frequency of least frequent allele
	allele_frequency = 0.0
	unique_kmers = 0
	missing_alleles = 0

	# determine the allele frequency (of the least frequent genotype allele)
	if 'AF' in info_fields:
		afs = [float(i) for i in info_fields['AF'].split(',')]
		# add frequency of reference allele
		frequencies = [1.0 - sum(afs)] + afs
		allele_frequency = max(min(frequencies),0.0)

	# determine minimum number of kmers that cover each allele
	if 'AK' in info_fields:
		aks = [int(i) for i in info_fields['AK'].split(',')]
		# determine minimum of unique kmers (exclude alleles not covered by paths)
		unique_kmers = min([i for i in aks if i >= 0])

	# determine number of missing alleles at variant position
	if 'MA' in info_fields:
		missing_alleles = int(info_fields['MA'])

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
				result[var_id].append(Variant(all_genotypes[(var_id, var_type)], likelihood, allele_frequency, unique_kmers, missing_alleles))
		else:
			for i in allele_ids:
				result[i].append(Variant(3, likelihood, 0.0, 0, 0))
	return result


class GenotypingStatistics:

	def __init__(self, quality, allele_frequency, unique_kmers, missing_alleles):

		# thresholds used to filter variant set
		self.quality = quality
		self.allele_frequency = allele_frequency
		self.unique_kmers = unique_kmers
		self.missing_alleles = missing_alleles

		# numbers of variants in baseline and callset
		self.total_baseline = 0
		self.total_intersection = 0

		# fractions of correct, wrong and untyped variants
		self.correct_all = 0
		self.wrong_all = 0
		self.not_typed_all = 0

		self.not_in_callset_all = 0
		self.absent_truth = 0

		# confusion matrix
		self.confusion_matrix = [[0 for x in range(4)] for y in range(4)]

	def print_stats_to_file(self, tsv_file):
		typed_all = max(self.correct_all + self.wrong_all, 1)
		assert self.total_baseline == self.correct_all + self.wrong_all + self.not_in_callset_all + self.not_typed_all

		non_ref_correct = (self.confusion_matrix[1][1] + self.confusion_matrix[2][2]) / max(1.0, float(sum(self.confusion_matrix[1][:-1]) + sum(self.confusion_matrix[2][:-1])))


		tsv_file.write('\t'.join([	str(self.quality), # quality
						str(self.allele_frequency), # allele freq
						str(self.unique_kmers), # unique kmer count
						str(self.missing_alleles), # nr of missing alleles
						str(self.total_baseline), # total baseline variants
						str(self.total_intersection), # intersection of callsets
						str(self.correct_all/float(typed_all)), # correct
						str(non_ref_correct), # correct non-ref
						str(self.wrong_all/float(typed_all)), # wrong
						str((self.not_typed_all+self.not_in_callset_all)/float(self.total_baseline if self.total_baseline is not 0 else 1)), # not typed
						str(self.correct_all),
						str(self.wrong_all),
						str(self.not_typed_all),
						str(self.not_in_callset_all),
						str(self.absent_truth)]) + '\n'
						)

	def print_matrix_to_file(self, txt_file):
		txt_file.write('\n###############################################################################################################################################\n')
		txt_file.write('#  Matrix for quality: '  + str(self.quality) +', allele frequency threshold: ' + str(self.allele_frequency) + ' , unique kmer count threshold: ' + str(self.unique_kmers) + ' and missing alleles threshold: ' + str(self.missing_alleles) + '  #\n')
		txt_file.write('  ###############################################################################################################################################\n\n')
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

	def nr_baseline_variants(self):
		return self._total_baseline

	def nr_callset_variants(self):
		return self._total_callset

	def get_samples(self):
		return self._samples

	def print_statistics(self, sample, varianttypes, quality, allele_freq, uk_count, missing_count, txt_files, tsv_files):
		"""
		Compute genotype concordance statistics for all variants with
		a quality at least quality and for which the allele frequency 
		of the least frequent allele is at least allele_freq.
		"""

		# make sure there is a txt and tsv file for each variant type
		assert len(txt_files) == len(tsv_files) == len(varianttypes)
		# keep variant statistics for all variant types
		statistics = { var_type : GenotypingStatistics(quality, allele_freq, uk_count, missing_count) for var_type in varianttypes }
		sample_id = self._samples.index(sample)
		# compare genotype predictions for all variant IDs
		for var_id, sample_to_genotype in self._baseline_variants.items():
			var_type = id_to_vartype(var_id)
			# ground truth genotype
			baseline_gt = sample_to_genotype[sample_id]
			if baseline_gt.get_binary_genotype() == 3:
				# if true genotype is unknown, skip position
				statistics[var_type].absent_truth += 1
				continue
			statistics[var_type].total_baseline += 1
			callset_predictions = self._callset_variants[var_id]
			if len(callset_predictions) > 0:
				# predicted genotype for the sample
				callset_gt = callset_predictions[sample_id]
				statistics[var_type].total_intersection += 1
				# if thresholds are not met, treat genotype as untyped
				if callset_gt.get_quality() < quality or callset_gt.get_allele_frequency() < allele_freq or callset_gt.get_unique_kmers() < uk_count or callset_gt.get_missing_alleles() < missing_count:
					print('untyped', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
					statistics[var_type].not_typed_all += 1
					continue
				statistics[var_type].confusion_matrix[baseline_gt.get_binary_genotype()][callset_gt.get_binary_genotype()] += 1
				if callset_gt.get_binary_genotype() == 3:
					# untyped
					statistics[var_type].not_typed_all += 1
					print('untyped', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
					continue
				if callset_gt.get_binary_genotype() == baseline_gt.get_binary_genotype():
					# correct genotype prediction
					statistics[var_type].correct_all += 1
					print('correct', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
				else:
					# wrong genotype prediction
					print('wrong', var_type.name, sample, var_id, callset_gt.get_binary_genotype(), baseline_gt.get_binary_genotype())
					statistics[var_type].wrong_all += 1
			else:
				# variant not in callset
				statistics[var_type].not_in_callset_all += 1
				statistics[var_type].confusion_matrix[baseline_gt.get_binary_genotype()][3] += 1
		for vartype in VariantType:
			assert (statistics[vartype].correct_all + statistics[vartype].wrong_all + statistics[vartype].not_typed_all == statistics[vartype].total_intersection)
#			assert (statistics[vartype].total_intersection + statistics[vartype].not_in_callset_all + statistics[vartype].absent_truth == statistics[vartype].total_baseline)
		# write statistics and confusion matrices to files
		for i,vartype in enumerate(varianttypes):
			statistics[vartype].print_stats_to_file(tsv_files[i])
			statistics[vartype].print_matrix_to_file(txt_files[i])

if __name__ == "__main__":

	# baseline: the file containing variants and true genotypes
	# callset: contains genotype predictions of genotyper	
	parser = argparse.ArgumentParser(prog='genotype-concordance-sample.py', description=__doc__)
	parser.add_argument('baseline', metavar='BASELINE', help='baseline VCF (ground truth).')
	parser.add_argument('callset', metavar='CALLSET', help='callset VCF (genotyped variants).')
	parser.add_argument('prefix', metavar='OUTFILE', help='prefix of the output file name.')
	parser.add_argument('samples', metavar='SAMPLES', help='comma separated list of samples to evaluate.')
	parser.add_argument('--qualities', default='0', metavar='GQ-THRESHOLDS', help='comma separated list of GQ-thresholds to consider (default: consider all variants regardless of quality).')
	parser.add_argument('--allele-frequencies', default='0', metavar='AF-THRESHOLDS', help='comma separated list of allele frequency thresholds to consider (default: consider all variants regardless of allele frequency.).')
	parser.add_argument('--unique-kmers', default='0', metavar='UK-THRESHOLDS', help='comma separated list of unique kmer counts for a variant (only works for VCFs with UK tag.).')
	parser.add_argument('--missing', default='0', metavar='MISSING-ALLELES_THRESHOLDS', help='comma separated list of missing allele counts for a variant (only works for VCFs with MA tag.).')
	parser.add_argument('--use-qual', default=False, action='store_true', help='use qualities in QUAL field instead of GQ fields (default: use GQ field).')
	args = parser.parse_args()

	# parse the thresholds to be considered
	quality_thresholds = [int(i) for i in args.qualities.split(',')]
	allele_freq_thresholds = [float(i) for i in args.allele_frequencies.split(',')]
	uk_thresholds = [int(i) for i in args.unique_kmers.split(',')]
	missing_thresholds = [int(i) for i in args.missing.split(',')]
	samples = args.samples.split(',')
	genotype_concordance = GenotypeConcordanceComputer(args.baseline, args.callset, samples, args.use_qual)

	
	# compute statistics for each sample
	for sample in samples:
		# prepare output files
		variantnames = []
		for vartype in VariantType:
			name = vartype.name.replace('_', '-')
			variantnames.append(name)
		tsv_files = [open(args.prefix + '_' + sample + '_' + name + '.tsv', 'w') for name in variantnames]
		txt_files = [open(args.prefix + '_' + sample + '_' + name + '.txt', 'w') for name in variantnames]
		# write headers to tsv files
		header = '\t'.join([	'quality',
					'allele_frequency',
					'unique_kmers',
					'missing_alleles',
					'total_baseline',
					'total_intersection',
					'correct_all',
					'correct_non-ref',
					'wrong_all',
					'not_typed_all',
					'nr_correct_all',
					'nr_wrong_all',
					'nr_not_typed_all',
					'nr_not_in_callset_all',
					'absent_truth'	]) + '\n'
		for tsv_file in tsv_files:
			tsv_file.write(header)
		vartypes = [vartype for vartype in VariantType]
		# compute statistics using all variants (regardless of thresholds)
		genotype_concordance.print_statistics(sample, vartypes, 0, 0.0, 0, 0, txt_files, tsv_files)
		# consider different thresholds on number of unique kmers
		for uk_count in uk_thresholds:
			for quality in quality_thresholds:
				if uk_count == quality == 0:
					continue
				genotype_concordance.print_statistics(sample, vartypes, quality, 0.0, uk_count, 0, txt_files, tsv_files)
		# consider different thresholds on allele frequencies
		for allele_freq in allele_freq_thresholds:
			if allele_freq == 0:
				continue
			genotype_concordance.print_statistics(sample, vartypes, 0, allele_freq, 0, 0, txt_files, tsv_files)
		# consider different thresholds on number of missing alleles
		for missing_count in missing_thresholds:
			if missing_count == 0:
				continue
			genotype_concordance.print_statistics(sample, vartypes, 0, 0.0, 0, missing_count, txt_files, tsv_files)
		# close files
		for i in range(len(tsv_files)):
			tsv_files[i].close()
			txt_files[i].close()
