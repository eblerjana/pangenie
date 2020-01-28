#!/usr/bin/python

import sys
from collections import defaultdict
from collections import namedtuple
import argparse
import math
from decimal import Decimal 
import vcf

Position = namedtuple("Position", "chrom position")
Genotype = namedtuple("Genotype", "alleles binary likelihood allele_freq unique_kmers")

def print_matrix(confusion_matrix):
	print("\t0\t1\t2\t./.")
	for true in range(4):
		line = str(true) + '\t'
		for geno in range(4):
			line += str(confusion_matrix[true][geno]) + '\t'
		print(line)

def extract_call(record, read_gl=False, read_qual=False):
	final_gt = None
	binary_gt = 3
	likelihood = 0
	allele_frequency = 0.0
	unique_kmers = 0
	if (record.samples[0]['GT'][0] != '.') and (record.samples[0]['GT'][-1] != '.'):
		genotype_string = record.samples[0]['GT'].replace('/', '|')
		allele_list = genotype_string.split('|')
		gt = (int(allele_list[0]), int(allele_list[1]))
		alleles = [record.REF] + record.ALT
		final_gt = set([str(alleles[gt[0]]), str(alleles[gt[1]])])

		# compute likelihood of genotype
		if read_gl:
			# check if GQ field is present
			if 'GQ' in record.FORMAT.split(':'):
				likelihood = int(record.samples[0]['GQ'])
		if read_qual:
			if record.QUAL is not None:
				likelihood = int(record.QUAL)

		# in case of bi-allelic SNP store binary gt
		if len(alleles) < 3:
			binary_gt = gt[0] + gt[1]
		else:
			binary_gt = -1

		# determine the allele frequency (of the least frequent genotype allele)
		if 'AF' in record.INFO:
			info_fields = record.INFO['AF']
			# add frequency of reference allele
			frequencies = [1.0 - sum(info_fields)] + info_fields
			allele_frequency = min(frequencies)

		# determine minimum number of kmers that cover each allele
		if 'AK' in record.INFO:
			info_fields = record.INFO['AK']
			# determine minimum of unique kmers (exclude alleles not covered by paths)
			unique_kmers = min([i for i in info_fields if i >= 0])

		# determine number of unique kmers
#		if 'UK' in record.FORMAT:
#			unique_kmers = int(record.samples[0]['UK'])

		

	return Position(record.CHROM, record.POS), Genotype(final_gt, binary_gt, likelihood, allele_frequency, unique_kmers)


def consider_position(record, snps, small_indels, midsize_indels, large_indels):
	alleles = [record.REF] + record.ALT
	varlen = max([len(a) for a in alleles])

	is_snp = record.is_snp
	is_small_indel = record.is_indel and varlen < 20
	is_large_indel = record.is_indel and varlen > 50
	is_midsize_indel = record.is_indel and (varlen >= 20) and (varlen <= 50)

	if snps and is_snp:
		return True
	if small_indels and is_small_indel:
		return True
	if large_indels and is_large_indel:
		return True
	if midsize_indels and is_midsize_indel:
		return True

	# if no flags set, consider all positions
	if not any([snps, small_indels, midsize_indels, large_indels]):
		return True
	else:
		return False

class GenotypeConcordanceComputer:
	def __init__(self, baseline_vcf, callset_vcf, use_qual, snps, small_indels, midsize_indels, large_indels):
		self._baseline_variants = {}
		self._callset_variants = defaultdict(list)
		self._duplicated_positions = defaultdict(lambda:False)
		self._total_baseline = 0
		self._total_callset = 0

		# read baseline variants
		for record in vcf.Reader(open(baseline_vcf, 'r')):
			if not consider_position(record, snps, small_indels, midsize_indels, large_indels):
				continue
			# extract position and genotype
			pos, genotype = extract_call(record)
			if pos in self._baseline_variants:
				# duplicated position, ignore
				print('Warning: position ' + str(pos.chrom) + ' ' + str(pos.position) + ' occurs more than once and will be skipped.')
				self._duplicated_positions[pos] = True
				continue
			self._total_baseline += 1
			self._baseline_variants[pos] = genotype

		# read callset variants
		for record in vcf.Reader(open(args.callset, 'r')):
			if not consider_position(record, snps, small_indels, midsize_indels, large_indels):
				continue
			if args.use_qual:
				pos, gt = extract_call(record, read_qual = True)
			else:
				pos, gt = extract_call(record, read_gl = True)
			self._callset_variants[pos].append(gt)
			self._total_callset += 1

	def print_statistics(self, quality, allele_freq, uk_count, tsv_output = None):
		"""
		Compute genotype concordance statistics for all variants with
		a quality at least quality and for which the allele frequency 
		of the least frequent allele is at least allele_freq.
		"""
		print('UK COUNT: ', uk_count)
		# confusion matrix for bi-allelic SNPs
		confusion_matrix = [[0 for x in range(4)] for y in range(4)]
		correct = 0
		wrong = 0
		no_prediction = 0
		not_in_callset = 0
		total_intersection = 0
		biallelic_snps = 0

		# check genotype predictions for each SNP
		for pos, gt in self._baseline_variants.items():

			# if position was present multiple times in baseline, skip it
			if self._duplicated_positions[pos]:
				continue

			callset_predictions = self._callset_variants[pos]

			if gt.binary != -1:
				biallelic_snps +=1 

			if len(callset_predictions) > 1:
				print('Warning: multiple predictions for variant at position ' + str(pos.chrom) + ' ' + str(pos.position) + '. Using first.')

			if len(callset_predictions) > 0:
				callset_gt = callset_predictions[0]
				total_intersection += 1

				if callset_gt.likelihood < quality or callset_gt.allele_freq < allele_freq or callset_gt.unique_kmers < uk_count:
					no_prediction += 1
					continue

				if (callset_gt.binary != -1) and (gt.binary!= -1):
					confusion_matrix[gt.binary][callset_gt.binary] += 1

				if (callset_gt.alleles == None) or (gt.alleles == None):
					no_prediction += 1
					continue

				# check if genotype predictions are identical
				if callset_gt.alleles == gt.alleles:
					correct += 1
				else:
					print('wrong', pos.chrom, pos.position, callset_gt.alleles, gt.alleles, gt.binary)
					wrong += 1
			else:
				not_in_callset += 1

		assert(correct + wrong + no_prediction == total_intersection)
		all_typed = correct + wrong if (correct + wrong) != 0 else 1
		no_prediction += not_in_callset
		total_non_duplicates_baseline = correct + wrong + no_prediction

		print('\n-----------------------------------------------------------------------------------------------------------------------------------------------------')
		print('                   STATISTICS for quality threshold: ' + str(quality) +', allele frequency threshold: ' + str(allele_freq) + ' and unique kmer count threshold: ' + str(uk_count))
		print('-----------------------------------------------------------------------------------------------------------------------------------------------------\n')
		print('matching genotypes (among all typed):\t' + str(correct/float(all_typed if all_typed != 0 else 1)) + ' ('+str(correct)+'/'+str(all_typed)+')')
		print('non-matching (among all typed):\t' + str(wrong/float(all_typed if all_typed != 0 else 1)) + ' ('+str(wrong)+'/'+str(all_typed)+')')
		print('variants genotyped as ./. in callset (among all variants not duplicated in baseline):\t' + str(no_prediction/float(total_non_duplicates_baseline)) + ' ('+str(no_prediction)+'/'+str(total_non_duplicates_baseline)+')')
		print('total ' + args.baseline + ':\t' + str(self._total_baseline))
		print('total ' + args.callset + ':\t' + str(self._total_callset))
		print('intersection:\t' + str(total_intersection))
		print('')
		print('total number of variants biallelic in baseline:\t' + str(biallelic_snps))
		print('CONFUSION MATRIX FOR BIALLELIC VARIANTS (vertical: ' + args.baseline + ' horizontal: ' + args.callset)
		print_matrix(confusion_matrix)
		typed_biallelic = 0
		for i in range(0,3):
			for j in range(0,3):
				typed_biallelic += confusion_matrix[i][j]
		correct_biallelic = confusion_matrix[0][0] + confusion_matrix[1][1] + confusion_matrix[2][2]
		wrong_biallelic = typed_biallelic - correct_biallelic
		no_prediction_biallelic = biallelic_snps - correct_biallelic - wrong_biallelic
		print('matching biallelic genotypes (among all typed):\t' + str(correct_biallelic/ float(typed_biallelic if typed_biallelic != 0 else 1)) + ' (' + str(correct_biallelic) + '/' + str(typed_biallelic) + ')')
		print('non-matching (among all typed):\t' + str(wrong_biallelic/ float(typed_biallelic if typed_biallelic != 0 else 1)) + ' (' + str(wrong_biallelic) + '/' + str(typed_biallelic) + ')')
		print('biallelic variants genotyped as ./. in callset (among all variants in intersection):\t' + str(no_prediction_biallelic/float(biallelic_snps)) + ' ('+str(no_prediction_biallelic)+'/'+str(biallelic_snps)+')')
		print('')

		if not tsv_output is None:
			tsv_output.write('\t'.join([str(quality), # quality
																	str(allele_freq), # allele freq
																	str(uk_count), # unique kmer count
																	str(correct/float(all_typed if all_typed != 0 else 1)), # correct
																	str(wrong/float(all_typed if all_typed != 0 else 1)), # wrong
																	str(no_prediction/float(total_non_duplicates_baseline)), # not typed
																	str(correct_biallelic/ float(typed_biallelic if typed_biallelic != 0 else 1)), # correct biallelic
																	str(wrong_biallelic/ float(typed_biallelic if typed_biallelic != 0 else 1)), # wrong biallelic
																	str(no_prediction_biallelic/float(biallelic_snps))  # not typed biallelic
														]) + '\n'
											) 


if __name__ == "__main__":

	# baseline: the file containing variants and true genotypes
	# callset: contains genotype predictions of genotyper	
	parser = argparse.ArgumentParser(prog='genotype-concordance.py', description=__doc__)
	parser.add_argument('baseline', metavar='BASELINE', help='baseline VCF (ground truth).')
	parser.add_argument('callset', metavar='CALLSET', help='callset VCF (genotyped variants).')
	parser.add_argument('--qualities', default='0', metavar='GQ-THRESHOLDS', help='comma separated list of GQ-thresholds to consider (default: consider all variants regardless of quality).')
	parser.add_argument('--allele-frequencies', default='0', metavar='AF-THRESHOLDS', help='comma separated list of allele frequency thresholds to consider (default: consider all variants regardless of allele frequency.).')
	parser.add_argument('--unique-kmers', default='0', metavar='UK-THRESHOLDS', help='comma separated list of unique kmer counts for a variant (only works for VCFs with UK tag.).')
	parser.add_argument('--use-qual', default=False, action='store_true', help='use qualities in QUAL field instead of GQ fields (default: use GQ field).')
	parser.add_argument('--snps', default=False, action='store_true', help='only consider SNPs.')
	parser.add_argument('--small-indels', default=False, action='store_true', help='only consider small indels ( length < 20bp).')
	parser.add_argument('--midsize-indels', default=False, action='store_true', help='only consider midsize indels ( 20 <= length <= 50bp )')
	parser.add_argument('--large-indels', default=False, action='store_true', help='only consider large indels ( length > 50 bp)')
	parser.add_argument('--tsv', metavar='TSV', default=None, help='store statistics in TSV-file.')
	args = parser.parse_args()

	quality_thresholds = [int(i) for i in args.qualities.split(',')]
	allele_freq_thresholds = [float(i) for i in args.allele_frequencies.split(',')]
	uk_thresholds = [int(i) for i in args.unique_kmers.split(',')]
	genotype_concordance = GenotypeConcordanceComputer(args.baseline, args.callset, args.use_qual, args.snps, args.small_indels, args.midsize_indels, args.large_indels)
	tsv_output = None
	if not args.tsv is None:
		tsv_output = open(args.tsv, 'w')
		tsv_output.write('\t'.join(['quality', 'allele_frequency', 'unique_kmers', 'correct_all', 'wrong_all', 'not_typed_all', 'correct_biallelic', 'wrong_biallelic', 'not_typed_biallelic']) + '\n')

	# print statistics for all quality and allele frequency thesholds
	genotype_concordance.print_statistics(0, 0.0, 0, tsv_output)

	print("\n######################################################################################################################################")
	print("                                          Evaluation by genotype quality thresholds")
	print("######################################################################################################################################")

	for uk_count in uk_thresholds:
		for quality in quality_thresholds:
			if uk_count == 0 and quality == 0:
				continue
			genotype_concordance.print_statistics(quality, 0.0, uk_count, tsv_output)

	print("\n######################################################################################################################################")
	print("                                          Evaluation by allele frequency thresholds")
	print("######################################################################################################################################")
	for allele_freq in allele_freq_thresholds:
		if allele_freq == 0:
			continue
		genotype_concordance.print_statistics(0, allele_freq, 0, tsv_output)



