#!/usr/bin/python

import argparse
import math
import sys
import random

parser = argparse.ArgumentParser(prog='sample-haplotypes.py', description=__doc__)
parser.add_argument('sample', metavar='sample', help='sample to genotype (will be excluded from panel, if included).')
parser.add_argument('number', metavar='nr_samples', help='number of samples to choose.', type=int)
parser.add_argument('panel', metavar='panel', help='name of output reference panel VCF file.')
parser.add_argument('--truth', metavar='truth', help='also generate ground truth VCF containing genotypes of given sample and write it to this file.', default=None)
args = parser.parse_args()

truth_file=None
if args.truth:
	# ground truth output file
	truth_file = open(args.truth, 'w')

selected_sample_ids=None
sample_id=None
panel_file=open(args.panel, 'w')

for line in sys.stdin:
	if line.startswith('##'):
		# description line
		if (args.truth):
			truth_file.write(line)
		panel_file.write(line)
		continue
	if line.startswith('#'):
		# header line
		if selected_sample_ids != None:
			raise ValueError('Malformatted input VCF file.')
		fields=line.split()
		# all vcf samples
		all_samples = fields[9:]
		# randomly choose args.number samples
		samples_to_choose=[]
		for i,s in enumerate(all_samples):
			if s == args.sample:
				sample_id = i + 9
			else:
				samples_to_choose.append( (s,i + 9) )
		if len(samples_to_choose) < args.number:
			raise ValueError('Not enough samples given.')
		selected_sample_ids=random.sample(samples_to_choose, args.number)
		selected_sample_ids.sort(key=lambda tup: tup[1])
		sys.stderr.write('Selected samples: ' + ','.join([s[0] for s in selected_sample_ids]) + '\n')
		# output header line
		header_line = fields[:9]
		panel_file.write('\t'.join(header_line + [s[0] for s in selected_sample_ids]) + '\n')
		if args.truth:
			truth_file.write('\t'.join(header_line) + '\t' + args.sample + '\n')
		continue
	# vcf line
	if selected_sample_ids == None:
		raise ValueError('Malformatted input VCF file.')
	fields=line.split()
	panel_line = fields[:9]
	for s in selected_sample_ids:
		panel_line.append(fields[s[1]])
	# if all genotypes are 0|0, skip line
	if all(x == '0|0' for x in panel_line[9:]):
		sys.stderr.write('Skipping position ' + fields[0] + ':' + fields[1] + ' since no alternative alleles are present in selected samples.\n')
		continue
	panel_file.write('\t'.join(panel_line) + '\n')
	# write truth file if requested
	if args.truth:
		if sample_id == None:
			raise RuntimeError('Truth set requested but given sample not contained in input panel.')
		truth_file.write('\t'.join(fields[:9] + [fields[sample_id]]) + '\n')	

if (args.truth):
	truth_file.close()
panel_file.close()
		
