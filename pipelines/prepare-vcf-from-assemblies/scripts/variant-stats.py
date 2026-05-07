import sys
from collections import defaultdict
from variantclassifier import VariantType, determine_variant_type

# output folder
output = sys.argv[1]

samples = []
stats = []

outfiles = []
chromosomes = ['chr' + str(i) for i in range(1,23)] + ['chrX']

numbers_varianttypes = defaultdict(int)
numbers_sample = defaultdict(int)
total_sample = defaultdict(int)

for line in sys.stdin:
	if line.startswith('##'):
		continue
	fields = line.strip().split()
	if line.startswith('#'):
		samples = fields[9:]
		stats = [0] * len(samples)
		outfiles = [open(output + '/' + sample + '-untypable.tsv', 'w') for sample in samples]
		continue
	if fields[0] not in chromosomes:
		continue
	present = []
	var_id = fields[2]
	var_type = determine_variant_type(line)
	assert var_id != '.'
	numbers_varianttypes[var_type] += 1
	for sample, genotype in zip(samples, fields[9:]):
		alleles = genotype.replace('|', '/').split('/')
		assert len(alleles) == 2
		assert alleles[0] in ['.', '0', '1']
		assert alleles[1] in ['.', '0', '1']
		if alleles[0] == "1" or alleles[1] == "1":
			present.append(sample)
			total_sample[(sample,var_type)] += 1
	assert len(present) > 0
	if len(present) == 1:
		numbers_sample[(present[0], var_type)] += 1
		sample_index = samples.index(present[0])
		stats[sample_index] += 1
		outfiles[sample_index].write(var_id + '\n')

for o in outfiles:
	o.close()

samples = sorted(set([s[0] for s in total_sample.keys()]))
variants = sorted([v for v in VariantType], key=lambda s: s.name)

header = ['type'] + [sample + '\t' for sample in samples] + ['total']
subheader = [''] + ['total','untypable']*len(samples) + ['\t']
print('\t'.join(header))
print('\t'.join(subheader))

for variant in variants:
	values = [variant.name]
	for sample in samples:
		values.append(str(total_sample[(sample, variant)]))
		values.append(str(numbers_sample[(sample, variant)]))
	values.append(str(numbers_varianttypes[variant]))
	print('\t'.join(values))

