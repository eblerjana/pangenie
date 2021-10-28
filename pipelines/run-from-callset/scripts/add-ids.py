import sys

allele_index = 0

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
	info_fields = { i.split('=')[0] : i.split('=')[1] for i in info_string.split(';') if '=' in i}
	alleles = fields[4].split(',')
	chrom = fields[0]
	position = fields[1]
	ids = []
	for i,a in enumerate(alleles):
		var_len = len(a)
		var_type = 'allele' + str(allele_index)
		var_id = '-'.join([chrom, position, var_type, str(var_len)])
		ids.append(var_id)
		allele_index += 1
	info_fields['ID'] = ','.join(ids)
	updated_info = []
	for k,v in info_fields.items():
		updated_info.append(k + '=' + v)
	fields[7] = ';'.join(updated_info)
	# if there is additional info present in sample field, remove it
	format = fields[8].split(':')
	if not 'GT' in format:
		raise Exception('Input VCF does not contain GT field at position ' + chrom + ':' + position + '.')
	fields[8] = 'GT'
	gt_index = format.index('GT')
	for i,sample in enumerate(fields[9:]):
		fields[9 + i] = sample.split(':')[gt_index]
		if not '|' in fields[ 9 + i]:
			raise Exception('Input VCF contains unphased genotype at position ' + chrom + ':' + position + '.')
	print('\t'.join(fields)) 
