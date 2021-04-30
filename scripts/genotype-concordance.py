
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
	def __init__(self, genotype, binary_genotype, variant_type : VariantType, is_nonref, quality=None, allele_frequency=None, unique_kmers=None, missing_alleles=None):
		"""
		Represents a genotyped variant.
		"""
		self._genotype=genotype
		self._binary_genotype=binary_genotype
		self._variant_type=variant_type
		self._quality=quality
		self._allele_frequency=allele_frequency
		self._unique_kmers=unique_kmers
		self._missing_alleles=missing_alleles
		self.nonref = is_nonref # true if the genotype is different from 0/0, i.e. the sample carries at least one non-reference allele

	def get_genotype(self):
		return self._genotype

	def get_binary_genotype(self):
		return self._binary_genotype

	def is_type(self, vartype : VariantType):
		if vartype == self._variant_type:
			return True
		else:
			return False

	def get_type(self):
		return self._variant_type

	def get_quality(self):
		return self._quality

	def get_allele_frequency(self):
		return self._allele_frequency

	def get_unique_kmers(self):
		return self._unique_kmers

	def get_missing_alleles(self):
		return self._missing_alleles

def determine_type_from_ids(ids):
	vartypes = [i.split('-')[2] for i in ids]
	if all([v == 'SNV' for v in vartypes]):
		return VariantType.snp

	# determine variant length
	allele_lengths = []
	for var_id in ids:
		if 'SNV' in var_id:
			continue
		assert var_id.split('-')[-2] in ['INS', 'DEL']
		length = int(var_id.split('-')[-1])
		allele_lengths.append(length)
	varlen = max(allele_lengths)

	is_insertion = 'INS' in vartypes and not 'DEL' in vartypes
	is_deletion = 'DEL' in vartypes and not 'INS' in vartypes

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex
	assert(False)

def determine_type_from_record(record):
	"""
	Determines the variant type.
	"""
	alleles = [record.REF] + record.ALT
	varlen = max([len(a) for a in alleles])

	if record.is_snp:
		return VariantType.snp

	is_deletion = record.var_subtype == 'del'
	is_insertion = record.var_subtype == 'ins'

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex
	assert(False)


def extract_call(record, samples, read_gl=False, read_qual=False):
	"""
	Extract genotype information from a VCF Record.
	"""
	
	# list containing Variant for each sample (in order of samples)
	result = [None] * len(samples)

	# allele frequency of least frequent allele
	allele_frequency = 0.0
	unique_kmers = 0
	missing_alleles = 0
	
	# determine the allele frequency (of the least frequent genotype allele)
	if 'AF' in record.INFO:
		info_fields = record.INFO['AF']
		# add frequency of reference allele
		frequencies = [1.0 - sum(info_fields)] + info_fields
		allele_frequency = max(min(frequencies),0.0)

	# determine minimum number of kmers that cover each allele
	if 'AK' in record.INFO:
		info_fields = record.INFO['AK']
		# determine minimum of unique kmers (exclude alleles not covered by paths)
		unique_kmers = min([i for i in info_fields if i >= 0])

	# determine number of missing alleles at variant position
	if 'MA' in record.INFO:
		missing_alleles = record.INFO['MA']
	
	# determine the type of variant
	variant_type = determine_type_from_record(record)
	
	for call in record.samples:
		# genotype alleles
		genotype_sequences = None
		# binary genotype
		binary_genotype = 3
		# genotype quality
		likelihood = 0
		is_nonref = True
		sample_name = call.sample
		if sample_name not in samples:
			continue
		genotype_string = call['GT'].replace('/', '|')
		if (genotype_string[0] != '.') and (genotype_string[-1] != '.'):
			genotype_list = genotype_string.split('|')
			genotype = (int(genotype_list[0]), int(genotype_list[1]))
			is_nonref = genotype[0] != 0 or genotype[1] != 0
			variant_alleles = [record.REF] + record.ALT
			genotype_sequences = set([str(variant_alleles[genotype[0]]), str(variant_alleles[genotype[1]])])

			# compute likelihood of genotype
			if read_gl:
				# check if GQ field is present
				if 'GQ' in record.FORMAT.split(':'):
					likelihood = int(call['GQ'])
			if read_qual:
				if record.QUAL is not None:
					likelihood = int(record.QUAL)

			# in case of bi-allelic SNP store binary gt
			if len(variant_alleles) < 3:
				binary_genotype = genotype[0] + genotype[1]
			else:
				binary_genotype = -1

		result[samples.index(sample_name)] = Variant(genotype_sequences, binary_genotype, variant_type, is_nonref, likelihood, allele_frequency, unique_kmers, missing_alleles)

	return Position(record.CHROM, record.POS), result


class GenotypingStatistics:

	def __init__(self, quality, allele_frequency, unique_kmers, missing_alleles):

		# thresholds used to filter variant set
		self.quality = quality
		self.allele_frequency = allele_frequency
		self.unique_kmers = unique_kmers
		self.missing_alleles = missing_alleles

		# numbers of variants in baseline and callset
		self.total_baseline = 0
		self.total_baseline_nonref = 0
		self.total_baseline_biallelic = 0
		self.total_intersection = 0

		# fractions of correct, wrong and untyped variants
		self.correct_all = 0
		self.correct_nonref = 0
		self.wrong_all = 0
		self.wrong_nonref = 0
		self.not_typed_all = 0
		self.not_typed_nonref = 0

		self.not_in_callset_all = 0
		self.not_in_callset_biallelic = 0
		self.not_in_callset_nonref = 0
		self.absent_truth = 0

		# confusion matrix
		self.confusion_matrix = [[0 for x in range(4)] for y in range(4)]

	def print_stats_to_file(self, tsv_file):
		typed_all = max(self.correct_all + self.wrong_all, 1)
		typed_nonref = max(self.correct_nonref + self.wrong_nonref, 1)
		assert self.total_baseline == self.correct_all + self.wrong_all + self.not_in_callset_all + self.not_typed_all
		assert self.total_baseline_nonref == self.correct_nonref + self.wrong_nonref + self.not_in_callset_nonref + self.not_typed_nonref
		correct_biallelic = self.confusion_matrix[0][0] + self.confusion_matrix[1][1] + self.confusion_matrix[2][2]
		wrong_biallelic = self.confusion_matrix[0][1] + self.confusion_matrix[0][2] + self.confusion_matrix[1][0] + self.confusion_matrix[1][2] + self.confusion_matrix[2][1] + self.confusion_matrix[2][0]
		not_typed_biallelic = self.total_baseline_biallelic - correct_biallelic - wrong_biallelic - self.not_in_callset_biallelic
		typed_biallelic = max(correct_biallelic + wrong_biallelic, 1)

		tsv_file.write('\t'.join([
			str(self.quality), # quality
			str(self.allele_frequency), # allele freq
			str(self.unique_kmers), # unique kmer count
			str(self.missing_alleles), # nr of missing alleles
			str(self.total_baseline), # total baseline variants
			str(self.total_baseline_biallelic), # total biallelic variants
			str(self.total_baseline_nonref), # total nonref variants
			str(self.total_intersection), # intersection of callsets
			str(self.correct_all/float(typed_all)), # correct
			str(self.wrong_all/float(typed_all)), # wrong
			str((self.not_typed_all+self.not_in_callset_all)/float(self.total_baseline if self.total_baseline != 0 else 1)), # not typed
			str(correct_biallelic/ float(typed_biallelic)), # correct biallelic
			str(wrong_biallelic/ float(typed_biallelic)), # wrong biallelic
			str((not_typed_biallelic+self.not_in_callset_biallelic)/float(self.total_baseline_biallelic if self.total_baseline_biallelic != 0 else 1)),  # not typed biallelic
			str(self.correct_nonref/float(typed_nonref)), # correct nonref
			str(self.wrong_nonref/float(typed_nonref)), # wrong nonref
			str((self.not_typed_nonref+self.not_in_callset_nonref)/float(self.total_baseline_nonref if self.total_baseline_nonref != 0 else 1)),  # not typed nonref
			str(self.correct_all),
			str(self.wrong_all),
			str(self.not_typed_all),
			str(self.not_in_callset_all),
			str(correct_biallelic),
			str(wrong_biallelic),
			str(not_typed_biallelic),
			str(self.not_in_callset_biallelic),
			str(self.correct_nonref),
			str(self.wrong_nonref),
			str(self.not_typed_nonref),
			str(self.not_in_callset_nonref)
			]) + '\n'
											)

	def print_matrix_to_file(self, txt_file):
		txt_file.write('\n###############################################################################################################################################\n')
		txt_file.write('#       Matrix for quality: '  + str(self.quality) +', allele frequency threshold: ' + str(self.allele_frequency) + ' , unique kmer count threshold: ' + str(self.unique_kmers) + ' and missing alleles threshold: ' + str(self.missing_alleles) + '#\n')
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
		self._baseline_variants = {}
		self._callset_variants = defaultdict(list)
		self._duplicated_positions = defaultdict(lambda:False)
		self._total_baseline = 0
		self._total_callset = 0
		self._samples = samples

		# read baseline variants
		for record in vcf.Reader(open(baseline_vcf, 'r')):
			# determine variant type
			pos, variants = extract_call(record, self._samples)
			if pos in self._baseline_variants:
				# duplicated position, ignore it
				sys.stderr.write('Warning: position ' + str(pos.chrom) + ' ' + str(pos.position) + ' occurs more than once and will be skipped.\n')
				self._duplicated_positions[pos] = True
				continue
			self._total_baseline += 1
			self._baseline_variants[pos] = variants

		# read callset variants
		for record in vcf.Reader(open(callset_vcf, 'r')):
			if use_qual:
				pos, variants = extract_call(record, self._samples, read_qual = True)
			else:
				pos, variants = extract_call(record, self._samples, read_gl = True)
			self._callset_variants[pos].append(variants)
			self._total_callset += 1

		sys.stderr.write('Read ' + str(self._total_baseline) + ' variants from the baseline VCF (including duplicates).\n')
		sys.stderr.write('Read ' + str(self._total_callset) + ' variants from the callset VCF (including duplicates).\n')

	def nr_baseline_variants(self):
		return self._total_baseline

	def nr_callset_variants(self):
		return self._total_callset

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
		# check genotype predictions
		for pos, genotypes in self._baseline_variants.items():

			gt = genotypes[sample_id]
			# if position was present multiple times in baseline, skip it
			if self._duplicated_positions[pos]:
				continue

			vartype = gt.get_type()

			# if true genotype is unknown, skip
			if gt.get_genotype() == None:
				statistics[vartype].absent_truth += 1
				continue

			statistics[vartype].total_baseline += 1

			if gt.nonref:
				statistics[vartype].total_baseline_nonref += 1

			callset_predictions = self._callset_variants[pos]

			if gt.get_binary_genotype() != -1:
				statistics[vartype].total_baseline_biallelic += 1

			if len(callset_predictions) > 1:
				print('Warning: multiple predictions for variant at position ' + str(pos.chrom) + ' ' + str(pos.position) + '. Using first.')

			if len(callset_predictions) > 0:
				callset_gt = callset_predictions[0][sample_id]
				statistics[vartype].total_intersection += 1

				if callset_gt.get_quality() < quality or callset_gt.get_allele_frequency() < allele_freq or callset_gt.get_unique_kmers() < uk_count or callset_gt.get_missing_alleles() < missing_count:
					statistics[vartype].not_typed_all += 1
					if gt.nonref:
						statistics[vartype].not_typed_nonref += 1
					continue

				if (callset_gt.get_binary_genotype() != -1) and (gt.get_binary_genotype() != -1):
					statistics[vartype].confusion_matrix[gt.get_binary_genotype()][callset_gt.get_binary_genotype()] += 1

				if (callset_gt.get_genotype() == None): # or (gt.get_genotype()==None):
					statistics[vartype].not_typed_all += 1
					if gt.nonref:
						statistics[vartype].not_typed_nonref += 1
					print(pos.position, vartype.name, callset_gt.get_genotype() == None, gt.get_genotype()==None)
					continue

				# check if genotype predictions are identical
				if callset_gt.get_genotype() == gt.get_genotype():
					statistics[vartype].correct_all += 1
					if gt.nonref:
						statistics[vartype].correct_nonref += 1
				else:
					print('wrong', vartype.name, pos.chrom, pos.position, sample, callset_gt.get_genotype(), gt.get_genotype(), gt.get_binary_genotype())
					statistics[vartype].wrong_all += 1
					if gt.nonref:
						statistics[vartype].wrong_nonref += 1
			else:
				statistics[vartype].not_in_callset_all += 1
				if gt.nonref:
					statistics[vartype].not_in_callset_nonref += 1
				if gt.get_binary_genotype() != -1:
					statistics[vartype].confusion_matrix[gt.get_binary_genotype()][3] += 1
					statistics[vartype].not_in_callset_biallelic += 1

		assert (statistics[vartype].correct_all + statistics[vartype].wrong_all + statistics[vartype].not_typed_all == statistics[vartype].total_intersection)
		# write statistics and confusion matrices to file
		for i,vartype in enumerate(varianttypes):
			statistics[vartype].print_stats_to_file(tsv_files[i])
			statistics[vartype].print_matrix_to_file(txt_files[i])

if __name__ == "__main__":

	# baseline: the file containing variants and true genotypes
	# callset: contains genotype predictions of genotyper	
	parser = argparse.ArgumentParser(prog='genotype-concordance.py', description=__doc__)
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
		header = '\t'.join([
				'quality',
				'allele_frequency',
				'unique_kmers',
				'missing_alleles',
				'total_baseline',
				'total_baseline_biallelic',
				'total_baseline_nonref',
				'total_intersection',
				'correct_all',
				'wrong_all',
				'not_typed_all',
				'correct_biallelic',
				'wrong_biallelic',
				'not_typed_biallelic',
				'correct_non-ref',
				'wrong_non-ref',
				'not_typed_non-ref',
				'nr_correct_all',
				'nr_wrong_all',
				'nr_not_typed_all',
				'nr_not_in_callset_all',
				'nr_correct_biallelic',
				'nr_wrong_biallelic',
				'nr_not_typed_biallelic',
				'nr_not_in_callset_biallelic',
				'nr_correct_non-ref',
				'nr_wrong_non-ref',
				'nr_not_typed_non-ref',
				'nr_not_in_callset_non-ref', 
	]) + '\n'
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
