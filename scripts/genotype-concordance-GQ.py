#!/usr/bin/python

import sys
from collections import defaultdict
import argparse
import math
from decimal import Decimal 
import vcf

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

	return (record.CHROM, record.POS), final_gt, binary_gt, likelihood

# baseline: the file containing variants and true genotypes
# callset: contains genotype predictions from forward-backward alg.		
parser = argparse.ArgumentParser(prog='genotype-concordance-GQ.py', description=__doc__)
parser.add_argument('baseline', metavar='BASELINE', help='baseline VCF (ground truth).')
parser.add_argument('callset', metavar='CALLSET', help='callset VCF (genotyped variants).')
parser.add_argument('--thresholds', default='0', metavar='GQ-THRESHOLDS', help='comma separated list of GQ-thresholds to consider (default: consider all variants regardless of quality).')
parser.add_argument('--use-qual', default=False, action='store_true', help='use qualities in QUAL field instead of GQ fields (default: use GQ field).')
parser.add_argument('--snps', default=False, action='store_true', help='only consider SNPs.')
parser.add_argument('--indels', default=False, action='store_true', help='only consider Indels.')
args = parser.parse_args()

all_SNPs = defaultdict(list)
total_baseline = 0
thresholds = [int(i) for i in args.thresholds.split(',')]
duplicated_positions = []

# read first callset
for record in vcf.Reader(open(args.baseline, 'r')):
	if args.snps:
		# check if record is a SNP
		if not record.is_snp:
			continue
	if args.indels:
		# check if record is an Indel
#		if not record.is_indel:
		if record.is_snp:
			continue
	# extract position and genotype
	pos, gt, binary, gl = extract_call(record)
	# if same position occurs multiple times, ignore all calls at this position
	if pos in all_SNPs:
		print('Warning: position ' + str(pos[0]) + ' ' + str(pos[1]) + ' occurs more than once and will be skipped.')
		duplicated_positions.append(pos)
		continue
	total_baseline += 1
	all_SNPs[pos] = [gt, binary]

all_callset_SNPs = []
# read second callset
for record in vcf.Reader(open(args.callset, 'r')):
	if args.snps:
		# check if record is a SNP
		if not record.is_snp:
			continue
	if args.indels:
		# check if record is an Indel
#		if not record.is_indel:
		if record.is_snp:
			continue
	if args.use_qual:
		pos, gt, binary, gl = extract_call(record, read_qual=True)
	else:
		pos, gt, binary, gl = extract_call(record, read_gl=True)
	all_callset_SNPs.append( (pos, gt, binary, gl) )

# compute statistics for each GQ-threshold
for threshold in thresholds:
	# confusion matrix for bi-allelic SNPs
	confusion_matrix = [[0 for x in range(4)] for y in range(4)]
	correct = 0
	wrong = 0
	no_prediction = 0
	total_callset = 0
	total_intersection = 0
	biallelic_snps = 0

	# check genotype predictions for each SNP
	for snp in all_callset_SNPs:
		pos, gt, binary, gl = snp
		total_callset += 1

		# if position was present multiple times in baseline, skip it
		if pos in duplicated_positions:
			continue

		# only consider calls contained in baseline
		if pos in all_SNPs.keys():
			total_intersection += 1

			if gl < threshold:
				no_prediction += 1
				continue

			if (binary != -1) and (all_SNPs[pos][1] != -1):
				biallelic_snps += 1
				confusion_matrix[all_SNPs[pos][1]][binary] += 1

			if (gt == None) or (all_SNPs[pos][0] == None):
				no_prediction += 1
				continue

			# check if genotype predictions are identical
			if gt == all_SNPs[pos][0]:
				#print('correct',pos[0],pos[1], gt,all_SNPs[pos])
				correct += 1
			else:
				print('wrong', pos[0],pos[1], gt, all_SNPs[pos][0], all_SNPs[pos][1])
				wrong += 1

	assert(correct + wrong + no_prediction == total_intersection)
	all_typed = correct + wrong #if (correct + wrong) != 0 else 1

	print('###### STATISTICS for GQ-threshold: ' + str(threshold) +' ######')
	print('matching genotypes (among all typed): ' + str(correct/float(all_typed if all_typed != 0 else 1)) + ' ('+str(correct)+'/'+str(all_typed)+')')
	print('non-matching (among all typed): ' + str(wrong/float(all_typed if all_typed != 0 else 1)) + ' ('+str(wrong)+'/'+str(all_typed)+')')
	print('variants genotyped as ./. in callset (among all variants in intersection): ' + str(no_prediction/float(total_intersection)) + ' ('+str(no_prediction)+'/'+str(total_intersection)+')')
	print('total ' + args.baseline + ': ' + str(total_baseline))
	print('total ' + args.callset + ': ' + str(total_callset))
	print('intersection: ' + str(total_intersection))
	print('')
	print('total number of variants biallelic in both sets: ' + str(biallelic_snps))
	print('CONFUSION MATRIX FOR BIALLELIC VARIANTS (vertical: ' + args.baseline + ' horizontal: ' + args.callset)
	print_matrix(confusion_matrix)
	print('')
	print('TODO: validate this script!')
