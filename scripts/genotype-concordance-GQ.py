#!/usr/bin/python

import sys
from collections import defaultdict
from collections import namedtuple
import argparse
import math
from decimal import Decimal 
import vcf

Position = namedtuple("Position", "chrom position")
Genotype = namedtuple("Genotype", "alleles binary likelihood")

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
	return Position(record.CHROM, record.POS), Genotype(final_gt, binary_gt, likelihood)

# baseline: the file containing variants and true genotypes
# callset: contains genotype predictions of genotyper	
parser = argparse.ArgumentParser(prog='genotype-concordance-GQ.py', description=__doc__)
parser.add_argument('baseline', metavar='BASELINE', help='baseline VCF (ground truth).')
parser.add_argument('callset', metavar='CALLSET', help='callset VCF (genotyped variants).')
parser.add_argument('--thresholds', default='0', metavar='GQ-THRESHOLDS', help='comma separated list of GQ-thresholds to consider (default: consider all variants regardless of quality).')
parser.add_argument('--use-qual', default=False, action='store_true', help='use qualities in QUAL field instead of GQ fields (default: use GQ field).')
parser.add_argument('--snps', default=False, action='store_true', help='only consider SNPs.')
parser.add_argument('--indels', default=False, action='store_true', help='only consider Indels.')
args = parser.parse_args()

baseline_variants = {}
callset_variants = defaultdict(list)
duplicated_positions = defaultdict(lambda:False)
total_baseline = 0
total_callset = 0
quality_thresholds = [int(i) for i in args.thresholds.split(',')]

# read baseline variants
for record in vcf.Reader(open(args.baseline, 'r')):
	if args.snps:
		if not record.is_snp:
			continue
	if args.indels:
		if not record.is_indel:
#		if record.is_snp:
			continue
	# extract position and genotype
	pos, genotype = extract_call(record)
	if pos in baseline_variants:
		# duplicated position, ignore
		print('Warning: position ' + str(pos.chrom) + ' ' + str(pos.position) + ' occurs more than once and will be skipped.')
		duplicated_positions[pos] = True
		continue
	total_baseline += 1
	baseline_variants[pos] = genotype

# read callset variants
for record in vcf.Reader(open(args.callset, 'r')):
	if args.snps:
		if not record.is_snp:
			continue
	if args.indels:
		if not record.is_indel:
#		if record.is_snp:
			continue
	if args.use_qual:
		pos, gt = extract_call(record, read_qual = True)
	else:
		pos, gt = extract_call(record, read_gl = True)
	callset_variants[pos].append(gt)
	total_callset += 1

# compute statistics for each GQ-threshold
for threshold in quality_thresholds:
	# confusion matrix for bi-allelic SNPs
	confusion_matrix = [[0 for x in range(4)] for y in range(4)]
	correct = 0
	wrong = 0
	no_prediction = 0
	not_in_callset = 0
	total_intersection = 0
	biallelic_snps = 0

	# check genotype predictions for each SNP
	for pos, gt in baseline_variants.items():

		# if position was present multiple times in baseline, skip it
		if duplicated_positions[pos]:
			continue

		callset_predictions = callset_variants[pos]

		if gt.binary != -1:
			biallelic_snps +=1 

		if len(callset_predictions) > 1:
			print('Warning: multiple predictions for variant at position ' + str(pos.chrom) + ' ' + str(pos.position) + '. Using first.')

		if len(callset_predictions) > 0:
			callset_gt = callset_predictions[0]
			total_intersection += 1

			if callset_gt.likelihood < threshold:
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
	all_typed = correct + wrong #if (correct + wrong) != 0 else 1
	no_prediction += not_in_callset
	total_non_duplicates_baseline = correct + wrong + no_prediction

	print('###### STATISTICS for GQ-threshold: ' + str(threshold) +' ######')
	print('matching genotypes (among all typed):\t' + str(correct/float(all_typed if all_typed != 0 else 1)) + ' ('+str(correct)+'/'+str(all_typed)+')')
	print('non-matching (among all typed):\t' + str(wrong/float(all_typed if all_typed != 0 else 1)) + ' ('+str(wrong)+'/'+str(all_typed)+')')
	print('variants genotyped as ./. in callset (among all variants not duplicated in baseline):\t' + str(no_prediction/float(total_non_duplicates_baseline)) + ' ('+str(no_prediction)+'/'+str(total_non_duplicates_baseline)+')')
	print('total ' + args.baseline + ':\t' + str(total_baseline))
	print('total ' + args.callset + ':\t' + str(total_callset))
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
