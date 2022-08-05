import sys

def determine_class_from_alleles(ref_allele, alt_allele):
	if len(ref_allele) == len(alt_allele) == 1:
		return 'SNV'
	if (len(ref_allele) > 1) and (len(alt_allele) == 1):
		return 'DEL'
	if (len(ref_allele) == 1) and (len(alt_allele) > 1):
		return 'INS'
	return 'COMPLEX'


def determine_allele_length(ref, allele, vartype):
	if vartype == 'SNV':
		return 1
	elif vartype == 'INS':
		return len(allele) - 1
	elif vartype == 'DEL':
		return len(ref) - 1
	elif vartype == 'COMPLEX':
		return len(allele)
	else:
		assert(False)

nr_variants = 0
nr_alleles = 0
counter = 0
for line in sys.stdin:
	if line.startswith('##'):
		print(line[:-1])
		continue
	if line.startswith('#'):
		print('##INFO=<ID=ID,Number=A,Type=String,Description="Variant IDs.">')
		print(line[:-1])
		continue
	fields = line.split()
	info_string = fields[7]
	info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in info_string.split(';') if '=' in i}
	alleles = fields[4].split(',')
	if len(alleles) > 1:
		raise Exception('Input VCF is not biallelic. Use bcftools norm -m-any to convert it into biallelic representation.')
	chrom = fields[0]
	vcf_id = fields[2] if fields[2] != '.' else None
	position = fields[1]
	ids = []
	for i,a in enumerate(alleles):
		vartype = determine_class_from_alleles(fields[3], a)
		var_len = determine_allele_length(fields[3], a, vartype)
		var_id = '-'.join([chrom, position + '_' + str(counter), vartype, str(var_len)])
		ids.append(var_id)
		counter += 1
		nr_alleles += 1
	info_fields['ID'] = ','.join(ids)
	updated_info = []
	for k,v in info_fields.items():
		updated_info.append(k + '=' + v)
	fields[7] = ';'.join(updated_info)
	# if there is additional info present in sample field, remove it
	format = fields[8].split(':')
	print('\t'.join(fields))
	nr_variants += 1
sys.stderr.write('Number of input variants: ' + str(nr_variants) + '\n')
sys.stderr.write('Number of input alleles: ' + str(nr_alleles) + '\n')
