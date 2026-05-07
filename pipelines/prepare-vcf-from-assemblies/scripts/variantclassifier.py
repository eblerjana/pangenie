from enum import Enum

class VariantType(Enum):
	snp = 0
	small_insertion = 1
	small_deletion = 2
	small_complex = 3
	midsize_insertion = 4
	midsize_deletion = 5
	midsize_complex = 6
	large_insertion = 7
	large_deletion = 8
	large_complex = 9
	
def determine_type_from_ids(ids):
	"""
	Determine type of each individual allele
	and return a list of VariantTypes.
	"""
	ids = set(ids)
	results = []
	for var in ids:
		var_type = var.split('-')[2]
		var_len = int(var.split('-')[-1])
		if var_type == "SNV":
			results.append(VariantType.snp)
			continue
		if var_len < 20:
			if var_type == "INS":
				results.append(VariantType.small_insertion)
			elif var_type == "DEL":
				results.append(VariantType.small_deletion)
			else:
				results.append(VariantType.small_complex)
			continue

		if var_len >= 20 and var_len < 50:
			if var_type == "INS":
				results.append(VariantType.midsize_insertion)
			elif var_type == "DEL":
				results.append(VariantType.midsize_deletion)
			else:
				results.append(VariantType.midsize_complex)
			continue

		if var_len >= 50:
			if var_type == "INS":
				results.append(VariantType.large_insertion)
			elif var_type == "DEL":
				results.append(VariantType.large_deletion)
			else:
				results.append(VariantType.large_complex)
	return results


def determine_variant_type(line):
	"""
	Determines the variant type based on
	the IDs of the variant alleles.
	"""
	fields = line.split()
	info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
	assert 'ID' in info_fields
	allele_ids = info_fields['ID'].split(',')
	# handle merged IDs
	all_ids = []
	for i in allele_ids:
		for j in i.split(':'):
			all_ids.append(j)
	all_ids = list(set(all_ids))
	assert len(all_ids) > 0
	
	if all(['SNV' in i for i in all_ids]) and all([ i.split('-')[1] == all_ids[0].split('-')[1]  for i in all_ids]):
		# all SNVs starting at the same base
		return VariantType.snp

	is_deletion = (len(all_ids) == 1) and 'DEL' in all_ids[0]
	is_insertion = (len(all_ids) == 1) and 'INS' in all_ids[0]

	alleles = [fields[3]] + [f for f in fields[4].split(',')]
	varlen = max([len(a) for a in alleles])
	if is_deletion or is_insertion:
		assert len(all_ids) == 1
		varlen = int(all_ids[0].split('-')[-1])

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex


def determine_pyvcf_type(record):
	"""
	Determine variant type from pyvcf
	record which does not have IDs.
	"""
	alleles = [record.REF] + record.ALT

	if record.is_snp:
		return VariantType.snp

	is_deletion = (len(record.ALT) == 1) and (len(record.REF) > 1) and (len(record.ALT[0]) == 1)
	is_insertion = (len(record.ALT) == 1) and (len(record.REF) == 1) and (len(record.ALT[0]) > 1)
	
	varlen = max([len(a) for a in alleles])
	if is_deletion or is_insertion:
		# insertions and deletions always have one base more
		varlen -= 1

	assert varlen >= 1

	if '<' in str(record.ALT[0]):
		is_deletion = str(record.ALT[0])[1:-1] == 'DEL'
		is_insertion = str(record.ALT[0])[1:-1] == 'INS'	

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex


def determine_variant_length(record):
	"""
	Determine variant length from
	pyvcf record.
	"""
	assert 'ID' in record.INFO
	allele_ids = record.INFO['ID']
	# handle merged IDs
	all_ids = []
	for i in record.INFO['ID']:
		for j in i.split(':'):
			all_ids.append(j)
	assert len(all_ids) > 0
	if len(all_ids) > 1:
		alleles = [record.REF] + record.ALT
		return max([len(a) for a in alleles])
	else:
		return int(all_ids[0].split('-')[-1])


def determine_variant_class(record):
	"""
	Determine variant class (SNV, INS, DEL)
	from pyvcf record.
	"""			
	assert 'ID' in record.INFO
	allele_ids = record.INFO['ID']
	# handle merged IDs
	all_ids = []
	for i in record.INFO['ID']:
		for j in i.split(':'):
			all_ids.append(j)
	assert len(all_ids) > 0
	
	if all(['SNV' in i for i in all_ids]) and all([ i.split('-')[1] == all_ids[0].split('-')[1]  for i in all_ids]):
		return 'SNV'
	
	is_deletion = (len(all_ids) == 1) and 'DEL' in all_ids[0]
	is_insertion = (len(all_ids) == 1) and 'INS' in all_ids[0]
	
	if is_deletion:
		return 'DEL'
	if is_insertion:
		return 'INS'
	return 'COMPLEX'


def determine_class_from_line(line):
	"""
	Determine variant class (SNV, INS, DEL)
	from VCF line.
	"""
	splitted = line.split()
	alleles = [splitted[3]] + [s for s in splitted[4].split(',')]

	if all([len(a) == 1 for a in alleles]):
		return 'SNV'

	if len(alleles) > 2:
		return 'COMPLEX'

	is_deletion = (len(alleles[0]) > 1) and (len(alleles[1]) == 1)
	is_insertion = (len(alleles[0]) == 1) and (len(alleles[1]) > 1)

	if is_deletion:
		return 'DEL'
	elif is_insertion:
		return 'INS'
	else:
		return 'COMPLEX'

def determine_variant_from_line(line):
	splitted = line.split()
	alleles = [splitted[3]] + [s for s in splitted[4].split(',')]

	if all([len(a) == 1 for a in alleles]):
		return VariantType.snp

	info_fields = { i.split('=')[0] : i.split('=')[1] for i in splitted[7].split(';') if "=" in i}
	assert 'ID' in info_fields
	ids = info_fields['ID']
	is_deletion = (len(alleles[0]) > 1) and (len(alleles[1]) == 1) and (len(alleles) < 3)
	is_insertion = (len(alleles[0]) == 1) and (len(alleles[1]) > 1) and (len(alleles) < 3)

	varlen = max([len(a) for a in alleles]) - 1

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex
