#!/usr/bin/python

import sys
from collections import defaultdict
from collections import namedtuple
import argparse
import math
from decimal import Decimal 
import vcf

class Variant:
	def __init__(self, chrom, start, ref, alt):
		self.chrom = chrom
		self.start = start
		self.ref = ref
		# order of alt alleles can be different
		alts = alt.split(',')
		alts.sort()
		self.alt = ','.join(alts)

	def __eq__(self, other):
		if self.chrom != other.chrom:
			return False
		if self.start != other.start:
			return False
		if self.ref != other.ref:
			return False
		if self.alt != other.alt:
			return False
		return True

	def __hash__(self):
		return hash((self.chrom, self.start, self.ref, self.alt))

	def __repr__(self):
		return "Variant({}, {}, {}, {})".format(self.chrom, self.start, self.ref, self.alt)

def parse_line(line):
	fields = line.split()
	chrom = fields[0]
	start = int(fields[1])
	ref = fields[3]
	alt = fields[4]
	sample_fields = { f:s for f,s in zip(fields[8].split(':'), fields[9].split(':')) }
	assert 'GT' in sample_fields
	not_typed = '.' in sample_fields['GT']
	return Variant(chrom, start, ref, alt), not not_typed

def remove_genotype(line):
	fields = line.split()
	format_field = fields[8].split(':')
	sample_field = fields[9].split(':')
	sample_field[format_field.index('GT')] = './.'
	fields[9] = ':'.join(sample_field)
	return '\t'.join(fields) + '\n'
	

parser = argparse.ArgumentParser(prog='compare-predictions.py', description=__doc__)
parser.add_argument('callset1', metavar='CALLSET1', help='first callset to compare to. Take all variant positions genotyped in this callset as baseline.')
parser.add_argument('callset2', metavar='CALLSET2', help='second callset.')
parser.add_argument('output', metavar='OUTPUT', help='second callset with all positions not typed in first callset set to ./.')
args = parser.parse_args()

# read first callset and store positions that are typed
baseline_positions = {}

for line in open(args.callset1, 'r'):
	if line.startswith('#'):
		continue
	variant, is_typed = parse_line(line)
	baseline_positions[variant] = is_typed

# read second callset and only keep genotypes that were called in first set
callset1_variants = len(baseline_positions)
callset2_variants = 0

not_in_callset = 0
unchanged = 0
genotype_removed = 0

with open(args.output, 'w') as outfile:
	for line in open(args.callset2, 'r'):
		if line.startswith('#'):
			outfile.write(line)
			continue
		variant, is_typed = parse_line(line)
		callset2_variants += 1
		if not variant in baseline_positions:
			print(variant + ' not in first callset.')
			not_in_callset += 1
			continue
		else:
			if baseline_positions[variant]:
				outfile.write(line)
				unchanged += 1
			else:
				# remove genotype information
				print('remove genotype information for variant ' + str(variant) + ' in output.')
				outfile.write(remove_genotype(line))
				genotype_removed += 1
print('total callset1: ' + str(callset1_variants))
print('total callset2: ' + str(callset2_variants))
print('variants in callset2 not in callset1:\t' + str(not_in_callset))
print('variants not changed in callset2:\t' + str(unchanged))
print('variants with genotypes removed:\t' + str(genotype_removed))
