import sys
import argparse

parser = argparse.ArgumentParser(prog='hwe.py', description=__doc__)
parser.add_argument('min_an', metavar='AN', help="minimum number of paths a variant must be covered.")
parser.add_argument('--chromosomes', metavar='CHROMOSOMES', nargs='+', default=[], help='Only select these chromosomes.')
args = parser.parse_args()

min_an = float(args.min_an)

total_variants = 0
total_alleles = 0
#total_lv = 0
#total_lv_alleles = 0
total_ac = 0
total_ac_alleles = 0
total_an = 0
total_an_alleles = 0
total_chrom = 0
total_chrom_alleles = 0
total_ns = 0
total_ns_alleles = 0
total_missing_alt = 0
total_missing_alt_alleles = 0
total_filtered = 0
total_filtered_alleles = 0

samples = []

for line in sys.stdin:
	fields = line.split()
	if line.startswith('##'):
		print(line[:-1])
		continue
	if line.startswith('#'):
		print('\t'.join(fields).strip())
		# keep samples
		samples = fields[9:]
		continue
	chromosome = fields[0].split('.')[-1]
	an_threshold = min_an * 2 * len(samples)
	select_chromosome = (chromosome in args.chromosomes) or (args.chromosomes == [])
	total_variants += 1
	info_fields = {a.split('=')[0] : a.split('=')[1].strip() for a in fields[7].split(';') if '=' in a}
	nr_alt_alleles = len(fields[4].split(','))
	total_alleles += nr_alt_alleles
	assert 'AN' in info_fields
	assert 'AC' in info_fields
	ac = sum([int(i) for i in info_fields['AC'].split(',')])
	an = int(info_fields['AN'])
	no_n = not 'N' in fields[3] and not 'N' in fields[4]
	no_missing_alt = all(c in 'CAGTcagt,' for c in fields[4])

	# replace "." genotype by ".|."
	for i in range(len(fields[9:])):
		if (fields[i+9] == '.'):
			fields[i+9] = '.|.'
	if no_n:
		total_ns += 1
		total_ns_alleles += nr_alt_alleles
	if ac > 0:
		total_ac += 1
		total_ac_alleles += nr_alt_alleles
	if an >= an_threshold:
		total_an += 1
		total_an_alleles += nr_alt_alleles
	if select_chromosome:
		total_chrom += 1
		total_chrom_alleles += nr_alt_alleles
	if no_missing_alt:
		total_missing_alt += 1
		total_missing_alt_alleles += nr_alt_alleles
	if (an >= an_threshold) and select_chromosome and no_n and no_missing_alt and (ac > 0):
		fields[0] = chromosome
		outline = '\t'.join(fields).strip()
		total_filtered += 1
		total_filtered_alleles += nr_alt_alleles
		print(outline)	

sys.stderr.write('total number of variants: ' + str(total_variants) + '\n')
sys.stderr.write('total number of alleles: ' + str(total_alleles)  + '\n')
sys.stderr.write('number of variants AN>=' + str(min_an) + ': ' + str(total_an)  + '\n')
sys.stderr.write('number of alleles AN>=' + str(min_an) +  ': ' + str(total_an_alleles)  + '\n')
sys.stderr.write('number of variants AC>0: ' + str(total_ac)  + '\n')
sys.stderr.write('number of alleles AC>0: ' + str(total_ac_alleles)  + '\n')
sys.stderr.write('number of variants on selected chromosomes: ' + str(total_chrom)  + '\n')
sys.stderr.write('number of alleles on selected chromosomes: ' + str(total_chrom_alleles)  + '\n')
sys.stderr.write('number of variants not containing Ns: ' + str(total_ns) + '\n')
sys.stderr.write('number of alleles not containing Ns: ' + str(total_ns_alleles) + '\n')
sys.stderr.write('number of variants not containing missing ALTs: ' + str(total_missing_alt) + '\n')
sys.stderr.write('number of alleles not containing missing ALTs: ' + str(total_missing_alt_alleles) + '\n')
sys.stderr.write('total number of variants in output (AN>=' + str(min_an) + ', on chromosomes): ' + str(total_filtered)  + '\n')
sys.stderr.write('total number of alleles in output: ' + str(total_filtered_alleles)  + '\n')
