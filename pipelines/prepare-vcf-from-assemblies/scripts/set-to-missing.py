import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(prog='set-to-missing.py', description="cat <vcf-file> | python3 set-to-missing.py ")
parser.add_argument('-v', metavar='VCF', required=True, help='The vcf file.')
parser.add_argument('-f', metavar='IDS', nargs='+', default=[], required=True, help='Variant IDs covered by the assemblies. One file per haplotype. Filename corresponds to the name of the sample.')
parser.add_argument('-m', metavar='MISSING', type=float, default=1.0, help='Skip positions if fraction of missing alleles is greater than this value.')
args = parser.parse_args()

# read the variant ids covered by each assembly
variant_to_sample = defaultdict(lambda: defaultdict(lambda: False))
# stats
sample_to_missing = defaultdict(lambda: 0)
total_variants = 0
skipped_variants = 0
skipped_missing = 0

for file in args.f:
	sample_name = file.split('/')[-1][:-4]
	count = 0
	for var_id in open(file, 'r'):
		variant_to_sample[var_id.strip()][sample_name] = True
		count += 1
	sys.stderr.write('Read ' + str(count) + ' IDs for sample ' + sample_name + '.\n')

index_to_sample = []

for line in open(args.v):
	if line.startswith('##'):
		print(line[:-1])
		continue
	if line.startswith('#'):
		# determine column index of each sample
		fields = line.split()
		for i,sample in enumerate(fields[9:]):
			index_to_sample.append(sample)
		print(line[:-1])
		continue
	fields = line.split()
	var_id = fields[2]
	alt_present = False
	nr_missing = 0
	for i,sample in zip(range(len(fields[9:])), index_to_sample):
		if not variant_to_sample[var_id][sample]:
			fields[9+i] = '.'
			sample_to_missing[sample] += 1
		if not fields[9+i] in ['.', '0']:
			alt_present = True
		if fields[9+i] == '.':
			nr_missing += 1
	total_variants += 1
	frac_missing = nr_missing / float(len(index_to_sample))
	assert frac_missing <= 1.0
	if frac_missing > args.m:
		skipped_missing += 1
		continue
	if not alt_present:
		skipped_variants += 1
		continue
	print('\t'.join(fields))
	
# print stats
for sample in index_to_sample:
	sys.stderr.write('Set ' + str(sample_to_missing[sample]) + '/' + str(total_variants) + ' alleles to missing for ' + sample + '.\n')
sys.stderr.write('Skipped ' + str(skipped_variants) + ' since no alternative allele was carried by the samples.\n')
sys.stderr.write('Skipped ' + str(skipped_missing) + ' since the fraction of missing alleles was too high.\n')
