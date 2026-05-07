import gzip
import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(prog='correct-sex-chromosomes.py', description=__doc__)
parser.add_argument('vcf', metavar='VCF', help='VCF file')
parser.add_argument('sample_info', metavar='SAMPLES', help="minimum number of paths by which a variant must be covered.")
args = parser.parse_args()


# read VCF once and count the number of non-'.' alleles on each haplotype
alleles_per_haplotype_chrX = defaultdict(lambda: [0,0])
alleles_per_haplotype_chrY = defaultdict(lambda: [0,0])
samples = []

processed = 0

for line in gzip.open(args.vcf, 'rt'):
	processed += 1
	if processed % 100000 == 0:
		sys.stderr.write('Processed ' + str(processed) + ' lines.\n')
	if line.startswith('##'):
		continue
	fields = line.strip().split()
	if line.startswith('#'):
		samples = fields[9:]
		continue
	if fields[0] not in ['chrX', 'X', 'chrY', 'Y']:
		continue
	for genotype, sample in zip(fields[9:], samples):
		if '/' in genotype:
			# ignore unphased genotypes
			continue
		alleles = genotype.split('|')
		assert len(alleles) < 3
		for i, allele in enumerate(alleles):
			if allele != '.':
				if fields[0] in ['chrX', 'X']:
					alleles_per_haplotype_chrX[sample][i] += 1
				else:
					alleles_per_haplotype_chrY[sample][i] += 1


# read the VCF again and correct genotypes on chrX and chrY 
sample_to_sex = {s.strip().split()[0] : s.strip().split()[1] for s in open(args.sample_info, 'r')}

sys.stderr.write('chrX\n')
for k,v in alleles_per_haplotype_chrX.items():
	if k not in sample_to_sex:
		continue
	sys.stderr.write(k + ' ' + str(v[0]) + ' ' + str(v[1]) + ' ' + sample_to_sex[k] + '\n')

sys.stderr.write('chrY\n')
for k,v in alleles_per_haplotype_chrY.items():
	if k not in sample_to_sex:
		continue
	sys.stderr.write(k + ' ' + str(v[0]) + ' ' + str(v[1]) + ' ' + sample_to_sex[k] + '\n')

for line in gzip.open(args.vcf, 'rt'):
	if line.startswith('##'):
		print(line.strip())
		continue
	fields = line.strip().split()
	if line.startswith('#'):
		samples = fields[9:]
		print(line.strip())
		continue
	if fields[0] not in ['chrX', 'X', 'chrY', 'Y']:
		print(line.strip())
		continue
	if fields[0] in ['chrX', 'X']:
		updated_genotypes = []
		for genotype, sample in zip(fields[9:], samples):
			if sample not in sample_to_sex:
				updated_genotypes.append(genotype)
				continue
			if sample_to_sex[sample] == '2':
				# female sample, leave genotype unchanged
				updated_genotypes.append(genotype)
			else:
				# male sample, duplicate the haplotype with most non-'.' alleles
				assert sample_to_sex[sample] == '1'
				assert (not '/' in genotype)
				index =  1 if (alleles_per_haplotype_chrX[sample][0] < alleles_per_haplotype_chrX[sample][1]) else 0
				new_alleles = [genotype.strip().split('|')[index], genotype.strip().split('|')[index]]
				updated_genotypes.append('|'.join(new_alleles))
		fields = fields[:9] + updated_genotypes
		print('\t'.join(fields))

	if fields[0] in ['chrY', 'Y']:
		updated_genotypes = []
		for genotype, sample in zip(fields[9:], samples):
			if sample not in sample_to_sex:
				updated_genotypes.append(genotype)
				continue
			if sample_to_sex[sample] == '2':
				# female sample, set genotype to 0|0 since females don't carry Y chromosomes
				updated_genotypes.append('0|0')
			else:
				# male sample, duplicate the haplotype with most non-'.' alleles
				assert sample_to_sex[sample] == '1'
				if '|' in genotype:
					index =  1 if (alleles_per_haplotype_chrY[sample][0] < alleles_per_haplotype_chrY[sample][1]) else 0
				else:
					assert '/' in genotype
					index = 0 if genotype.split('/')[0] != '.' else 1
				new_alleles = [genotype.strip().replace('/', '|').split('|')[index], genotype.strip().replace('/', '|').split('|')[index]]
				updated_genotypes.append('|'.join(new_alleles))
		fields = fields[:9] + updated_genotypes
		print('\t'.join(fields))

