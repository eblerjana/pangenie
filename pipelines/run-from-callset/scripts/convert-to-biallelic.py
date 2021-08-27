#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(prog='convert-to-biallelic.py', description='cat <multiallelic VCF> | python convert-to-biallelic.py <biallelic VCF>')
parser.add_argument('vcf', metavar='VCF', help='original VCF containing REF/ALT of each Variant ID.')
args = parser.parse_args()

# chromosome ->  ID -> [start, REF, ALT] per chromosome
chrom_to_variants = defaultdict(lambda: defaultdict(list))

# read the biallelic VCF containing REF/ALT for all variant IDs and store them
for line in open(args.vcf, 'r'):
	if line.startswith('#'):
		continue
	fields = line.split()
	info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';')}
	assert 'ID' in info_field
	ids = info_field['ID'].split(',')
	alleles = fields[4].split(',')
	for id,allele in zip(ids, alleles):
		chrom_to_variants[fields[0]][id] = [fields[1], fields[3], allele]

for line in sys.stdin:
	if line.startswith('#'):
		# header line
		if any([i in line for i in ['INFO=<ID=AF', 'INFO=<ID=AK', 'FORMAT=<ID=GL', 'FORMAT=<ID=KC']]):
			# these fields will not be contained in biallelic VCF
			continue
		print(line[:-1])
		continue
	fields = line.split()
	assert len(fields) > 7
	# parse the INFO field
	info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';')}
	assert 'ID' in info_field
	# determine ID string belonging to each allele (keep empty string for REF, as it does not have an ID)
	allele_to_ids = [''] + info_field['ID'].split(',')
	# collect all variant IDs in this region
	ids = set([])
	for i in info_field['ID'].split(','):
		for j in i.split(':'):
			ids.add((j,int(chrom_to_variants[fields[0]][j][0])))
	# sort the ids by the starting coordinate (to ensure the VCF is sorted)
	ids = list(ids)
	ids.sort(key=lambda x : x[1])
	# create a single, biallelic VCF record for each ID
	for (var_id, coord) in ids:
		vcf_line = fields[:9]
		# set start coordinate
		vcf_line[1] = str(coord)
		# set REF
		vcf_line[3] = chrom_to_variants[fields[0]][var_id][1]
		# set ALT
		vcf_line[4] = chrom_to_variants[fields[0]][var_id][2]
		# set INFO
		vcf_line[7] = 'ID=' + var_id
		# also add other INFO fields (except ID which was replaced)
		for k,v in info_field.items():
			if k == 'ID':
				continue
			if k in ['MA', 'UK']:
				values = ';' + k + '=' + v
				vcf_line[7] = vcf_line[7] + values
		# keep only GT and GQ
		vcf_line[8] = 'GT:GQ'
		# determine the genotype of each sample
		for sample_field in fields[9:]:
			# determine position of GT and GQ from FORMAT
			assert 'GT' in fields[8]
			assert 'GQ' in fields[8]
			format_field = fields[8].split(':')
			index_of_gt = format_field.index('GT')
			index_of_gq = format_field.index('GQ')
			genotype = sample_field.split(':')
			biallelic_genotype = []
			for allele in genotype[index_of_gt].replace('|', '/').split('/'):
				if allele == '.':
					# missing allele
					biallelic_genotype.append('.')
				else:
					if var_id in allele_to_ids[int(allele)].split(':'):
						biallelic_genotype.append('1')
					else:
						biallelic_genotype.append('0')
			vcf_line.append('/'.join(biallelic_genotype) + ':' + genotype[index_of_gq])
		print('\t'.join(vcf_line))
