#!/usr/bin/python

import argparse
from cyvcf2 import VCF
import sys
from collections import namedtuple
from collections import defaultdict

AlleleStats = namedtuple('AlleleStats','af ac an untyped')
GenotypeStats = namedtuple('GenotypeStats', 'heterozygosity het total')

def compute_allele_statistics(record):
	"""
	Compute allele related statistics.
	"""
	an = 0
	ac = 0
	unknown = 0
	for genotype in record.genotypes:
		alleles = genotype[:-1]
		assert 1 <= len(alleles) <= 2
		if len(alleles) == 1:
			# haploid genotype
			alleles.append(alleles[0])
		for a in alleles:
			if a == -1:
				unknown += 1
				continue
			assert a in [0,1]
			an += 1
			ac += a
	if an < 1:
		assert ac < 1
	af = ac / max(1.0, float(an))
	return AlleleStats(str(af), str(ac), str(an), str(unknown))


def read_uk(record):
	return str(record.INFO['UK'])


def compute_genotype_statistics(record, qualities=None):
	"""
	Compute genotype related statistics.
	"""
	counts = defaultdict(int)
	het_genotypes = 0
	total_genotypes = 0
	gqs = record.format('GQ') if qualities is not None else [None]*len(record.genotypes)
	for genotype, quality in zip(record.genotypes, gqs):
		alleles = genotype[:-1]
		assert 1 <= len(alleles) <= 2
		if len(alleles) == 1:
			# haploid genotype
			alleles.append(alleles[0])
		if not -1 in alleles:
			total_genotypes += 1
			if sum(alleles) == 1:
				assert 0 in alleles
				assert 1 in alleles
				het_genotypes += 1
			# read GQ
			if qualities is not None:
				for q in qualities:
					if int(quality) >= q:
						counts[q] += 1
	genotype_stats = GenotypeStats( str(het_genotypes / max(1.0, float(total_genotypes))), str(het_genotypes), str(total_genotypes))
	return genotype_stats, counts


parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description=__doc__)
parser.add_argument('panel', metavar='panel', help='biallelic panel variants.')
parser.add_argument('pangenie', metavar='pangenie', help='PanGenie biallelic genotyped variants for all unrelated samples.')
parser.add_argument('pangenie_all', metavar='pangenie_all', help='PanGenie biallelic genotyped variants for all samples (including children).')
args = parser.parse_args()

# compute statistics for input panel VCF
panel_reader = VCF(args.panel)
panel_samples = panel_reader.samples
panel_stats = {}

for variant in panel_reader:
	# require bi-allelic vcf with IDs
	assert len(variant.ALT) == 1
	var_id = variant.INFO['ID']
	allele_stats = compute_allele_statistics(variant)
	panel_stats[var_id] = allele_stats

sys.stderr.write('Done with panel.\n')


# compute statistics for all unrelated samples
pangenie_reader = VCF(args.pangenie)
pangenie_samples = pangenie_reader.samples
pangenie_stats = {}
quals = [0,200]

for variant in pangenie_reader:
	assert len(variant.ALT) == 1
	var_id = variant.INFO['ID']
	allele_stats = compute_allele_statistics(variant)
	genotype_stats, counts = compute_genotype_statistics(variant, quals)
	uk = read_uk(variant)
	pangenie_stats[var_id] = [allele_stats, genotype_stats, uk, counts]

sys.stderr.write('Done with unrelated.\n')

# compute statistics for all samples (unrelated + related)
pangenie_all_reader = VCF(args.pangenie_all)
pangenie_all_samples = pangenie_all_reader.samples
pangenie_all_stats = {}

for variant in pangenie_all_reader:
	assert len(variant.ALT) == 1
	var_id = variant.INFO['ID']
	allele_stats = compute_allele_statistics(variant)
	genotype_stats, counts = compute_genotype_statistics(variant, quals)
	uk = read_uk(variant)
	pangenie_all_stats[var_id] = [allele_stats, genotype_stats, uk, counts]


# print stats for all IDs in genotypes VCF
header = [ 	'variant_id',
		'panel_allele_freq',
		'panel_alternative_alleles',
		'panel_total_alleles',
		'panel_unknown_alleles',

		'pangenie-unrelated_allele_freq',
		'pangenie-unrelated_alternative_alleles',
		'pangenie-unrelated_total_alleles',
		'pangenie-unrelated_unknown_alleles',
		'pangenie-unrelated_heterozygosity',
		'pangenie-unrelated_heterozygous_genotypes',
		'pangenie-unrelated_total_genotypes',
		'pangenie-unrelated_unique_kmers',

		'pangenie-all_allele_freq',
		'pangenie-all_alternative_alleles',
		'pangenie-all_total_alleles',
		'pangenie-all_unknown_alleles',
		'pangenie-all_heterozygosity',
		'pangenie-all_heterozygous_genotypes',
		'pangenie-all_total_genotypes',
		'pangenie-all_unique_kmers'
	]

for q in quals:
	header.append('pangenie-unrelated_GQ>=' + str(q))
	header.append('pangenie-all_GQ>=' + str(q))

print('\t'.join(header))

assert len(pangenie_stats) == len(pangenie_all_stats)

for var_id in pangenie_stats:
	if not var_id in panel_stats:
		continue
		
	line = [	var_id,
			panel_stats[var_id].af,
			panel_stats[var_id].ac,
			panel_stats[var_id].an,
			panel_stats[var_id].untyped,

			pangenie_stats[var_id][0].af,
			pangenie_stats[var_id][0].ac,
			pangenie_stats[var_id][0].an,
			pangenie_stats[var_id][0].untyped,
			pangenie_stats[var_id][1].heterozygosity,
			pangenie_stats[var_id][1].het,
			pangenie_stats[var_id][1].total,
			pangenie_stats[var_id][2],

			pangenie_all_stats[var_id][0].af,
			pangenie_all_stats[var_id][0].ac,
			pangenie_all_stats[var_id][0].an,
			pangenie_all_stats[var_id][0].untyped,
			pangenie_all_stats[var_id][1].heterozygosity,
			pangenie_all_stats[var_id][1].het,
			pangenie_all_stats[var_id][1].total,
			pangenie_all_stats[var_id][2]

		]
		
	# add counts for GQs
	for q in quals:
		line.append(str(pangenie_stats[var_id][3][q]))
		line.append(str(pangenie_all_stats[var_id][3][q]))
	assert len(line) == len(header)
	print('\t'.join(line))
