import argparse
import gzip
import re
from collections import defaultdict


class Node:
	"""
	Represents a node in the GFA graph (node ID, direction).
	"""
	def __init__(self, id, direction):
		self.id = str(id)
		self.direction = direction
		self._marked = False
	def __hash__(self):
		return hash((self.id, self.direction))
	def __str__(self):
		return "(" + self.id + "," + self.direction + ')'
	def __eq__(self, other):
		if self.id != other.id:
			return False
		if self.direction != other.direction:
			return False
		return True
	def mark_node(self):
		self._marked = True
	def is_marked(self):
		return self._marked


class ReferencePath:
	"""
	Represents a reference traversal.
	"""
	def __init__(self, traversal):
		self._traversal = traversal
		self._start_node = traversal[0]
		self._end_node = traversal[-1]
		self._edges = {}
		for i in range(1, len(traversal)):
			self._edges[traversal[i-1]] = traversal[i]
	def contains_edge(self, start, end):
		if not start in self._edges:
			return False
		if self._edges[start] == end:
			return True
		else:
			return False
	def contains_node(self, node):
		if node == self._start_node or node == self._end_node:
			return True
		if node in self._edges:
			return True
		return False
	def __str__(self):
		return ",".join([str(a) for a in self._traversal])
	def get_subpath(self, start, end):
		subpath = []
		assert self.contains_node(start)
		assert self.contains_node(end)
		subpath.append(start)
		current = start
		while current != end:
			current = self._edges[current]
			subpath.append(current)
		return subpath
	def __eq__(self, other):
		if self._traversal != other._traversal:
			return False
		return True


def parse_gfa(filename, segments = None):
	"""
	Read nodes and their sequences from GFA
	and store them.
	"""
	node_to_data = {}
	with gzip.open(filename, 'rt') as gfa_file:
		for line in gfa_file:
			if not line[0] in ['S']:
				# we are only interested in the segments
				continue
			fields = line.split()
			if segments and fields[1] not in segments:
				continue
			ref_pos = None

			# search for SO:i entry
			so_entry = None
			for f in fields:
				if f.startswith('SO:i'):
					so_entry = f
			if len(fields) > 4 and so_entry:
				ref_pos = so_entry.split(':')[-1]
			node_to_data[fields[1]] = (ref_pos, fields[2])
	return node_to_data



def define_id(ref_allele, alt_allele, chrom, position, index):
	"""
	Define an unique identifier for a variant allele.
	"""
	len_ref = len(ref_allele)
	len_alt = len(alt_allele)
	vartype = None
	varlen = None
	if len_ref == 1:
		if len_alt > 1:
			vartype = 'INS'
			varlen = str(abs(len_ref - len_alt))
		else:
			vartype = 'SNV'
			varlen = '1'
	else:
		if len_alt == 1:
			vartype = 'DEL'
			varlen = str(abs(len_ref - len_alt))
		else:
			vartype = 'COMPLEX'
			varlen = str(max(len_ref, len_alt))
	return '-'.join([chrom, str(position), vartype, str(index), varlen])


def detect_variants(ref_traversal, alt_traversal):
	"""
	Given the reference traversal and an ALT traversal,
	determine the variants contained in the ALT traversal
	relative to the reference.
	"""
	node_to_index = defaultdict(list)
	for i,node in enumerate(alt_traversal):
		node_to_index[node].append(i)

	# if first and last nodes don't match, keep
	# allele and don't decompose
	if ref_traversal[0] != alt_traversal[0]:
		return [alt_traversal], False
	if ref_traversal[-1] != alt_traversal[-1]:
		return [alt_traversal], False

	prev_alt_index = 0
	prev_ref_index = 0
	alleles = []
	alt_end = len(alt_traversal) - 1
	ref_end = len(ref_traversal) - 1
	for i,node in enumerate(ref_traversal):
		if node in node_to_index:
			# find next occurance of the node
			alt_index = prev_alt_index
			for j in node_to_index[node]:
				if j >= alt_index:
					alt_index = j
					break
			if alt_index == prev_alt_index:
				# node is there, but is part of an earlier ALT allele (back edge)
				continue
			if (i == ref_end) and (alt_index != alt_end):
				# we hit the last reference node and that should always me
				# matched with the last node of the alt allele
				alt_index = alt_end
			if (abs(alt_index - prev_alt_index) > 1) or (abs(i-prev_ref_index) > 1):
				allele = alt_traversal[prev_alt_index : alt_index+1]
				if len(allele) > 1:
					alleles.append(allele)
			prev_alt_index = max(alt_index, prev_alt_index)
			prev_ref_index = i
	#validate(ref_traversal, alleles, alt_traversal)
	return alleles, True

# TODO: this is not a proper way to evaluate!
def validate(ref_traversal, alleles, alt_traversal):
	"""
	Given an allele traversal and its nested alleles relative
	to the reference traversal, check if inserting alleles in reference
	result in the alt_traversal.
	"""
	ref_string = traversal_to_string(ref_traversal)
	ref_path = ReferencePath(ref_traversal)
	alt_string = traversal_to_string(alt_traversal)

	expected_string = ref_string

	for allele in alleles:
		allele_string = traversal_to_string(allele)
		ref_str = traversal_to_string(ref_path.get_subpath(allele[0], allele[-1]))
		if ref_str not in expected_string:
			print(ref_string)
			print(allele_string)
		assert ref_str in expected_string
		assert expected_string.count(ref_str) == 1
		expected_string = expected_string.replace(ref_str, allele_string)
	if expected_string != alt_string:
		print(ref_string)
		print("ALLELE")
		print(alt_string)
		print("COMPUTED")
		print(expected_string)
		for allele in alleles:
			print(traversal_to_string(allele))
	assert expected_string == alt_string


def parse_allele_traversal(traversal):
	"""
	Converts a traversal string to a list of nodes.
	"""
	ids = re.split('<|>', traversal.strip())
	directions = ['>' if i == '>' else '<' for i in traversal if i in ['>', '<']]
	assert ids[0] == ''
	assert len(ids[1:]) == len(directions)
	allele_traversal = []
	for node, direction in zip(ids[1:], directions):
		allele_traversal.append(Node(node, direction))
	return allele_traversal


def traversal_to_string(traversal):
	"""
	Converts list of nodes representing a traversal
	back to a traversal string.
	"""
	result = ""
	for node in traversal:
		result += node.direction
		result += node.id
	return result



def parse_info(fields):
	"""
	Parse the info field of a VCF file
	"""
	info_fields = { k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k }
	return info_fields


def info_to_string(info_fields):
	"""
	Convert map containing INFO fields back to string
	"""
	return ';'.join([k + '=' + v for k,v in info_fields.items()])


def reverse_complement(sequence):
	"""
	Computes the reverse complement of a sequence.
	"""
	result = ''
	complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n'}
	i = len(sequence) - 1
	while i >= 0:
		result += complement[sequence[i]]
		i -= 1
	return result



def construct_allele_string(traversal, gfa, add_flank=True):
	"""
	Constructs the nucleotide sequence defined by
	an traversal.
	"""
	result = ''
	# first and last node are the same in any allele. Only take last base
	# of first one and skip the last one
	for i,node in enumerate(traversal[:-1]):
		assert node.id in gfa
		sequence = gfa[node.id][1] if node.direction == '>' else reverse_complement(gfa[node.id][1])
		if (i == 0):
			result += sequence[-1] if add_flank else ''
		else:
			result += sequence
	return result


def get_ref_position(node, gfa, add_flank=True):
	"""
	Look up reference position of a segment in graph.
	"""
	assert node.id in gfa
	position = int(gfa[node.id][0]) + len(gfa[node.id][1])
	if add_flank:
		return str(position)
	else:
		return str(position +1)


def trim_alleles(ref_pos, ref, alt):
	pos = int(ref_pos)
	while len(ref) >= 2 and len(alt) >= 2 and ref[-1] == alt[-1]:
		ref, alt = ref[:-1], alt[:-1]

	while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
		ref, alt = ref[1:], alt[1:]
		pos += 1
	return str(pos), ref, alt


def decompose(line, gfa):
	"""
	decomposes a large bubble into the ones nested inside (if any).
	"""
	fields = line.split()
	info_fields = parse_info(fields)
	assert 'AT' in info_fields
#	assert 'LV' in info_fields
#	if int(info_fields['LV']) > 0:
#		return None, None
	allele_traversals = info_fields['AT'].split(',')
	nr_alleles = len(allele_traversals)
	assert nr_alleles > 1
	# if a biallelic record assign ID and print directly
	if nr_alleles == 2:
		# biallelic variant only assign ID
		new_id = define_id(fields[3], fields[4], fields[0], fields[1], allele_traversals[1])
		info_fields['ID'] = new_id
		updated_info = info_to_string(info_fields)
		fields[7] = updated_info
		updated_line = '\t'.join(fields)
		return updated_line, [updated_line]
	else:
		ref_traversal = parse_allele_traversal(allele_traversals[0])
		ref_path = ReferencePath(ref_traversal)
		id_to_index = defaultdict(list)
		allele_to_ids = defaultdict(list)
		id_to_alleles = {}
		biallelic_records = []
		seen_variants = defaultdict(str)
		vcf_alleles = fields[4].split(',')

		# deconstruct alt alleles
		for i,a in enumerate(allele_traversals[1:]):
			alt_traversal = parse_allele_traversal(a)
			alleles, is_decomposed = detect_variants(ref_traversal, alt_traversal)
			if not is_decomposed:
				print('Could not decompose allele ' + str(i) + ' at position ' + fields[0] + ':' + fields[1])
				assert len(alleles) == 1
			for allele in alleles:
				if is_decomposed:
					ref_allele = ref_path.get_subpath(allele[0], allele[-1])
					add_flank = (len(allele) == 2) or (len(ref_allele) == 2)
					# translate traversal to string based on sequence information in graph
					alt_string = construct_allele_string(allele, gfa, add_flank)
					# determine reference allele and its sequence/position
					ref_string = construct_allele_string(ref_allele, gfa, add_flank)
					# determine reference position
					ref_pos = get_ref_position(allele[0], gfa, add_flank)
				else:
					ref_allele = ref_traversal
					ref_string = fields[3]
					alt_string = vcf_alleles[i]
					ref_pos = fields[1]
				# trim alleles (i.e. remove common prefix/suffix)
				ref_pos, ref_string, alt_string = trim_alleles(ref_pos, ref_string, alt_string)
				# in case allele strings are the same, there is no variant
				if ref_string == alt_string:
					print('Same allele sequence observed for ALT/REF traversals ' + traversal_to_string(allele) + ' and ' + traversal_to_string(ref_allele))
					continue
				if (ref_string, alt_string, ref_pos) in seen_variants:
					# if same alt allele was seen before, look up its allele_id
					allele_id = seen_variants[(ref_string, alt_string, ref_pos)]
					other_id = define_id(ref_string, alt_string, fields[0], ref_pos, traversal_to_string(allele))
					if other_id not in id_to_index:
						print('Sequence observed for ' + traversal_to_string(allele) + ' is the same as for ' + allele_id)
				else:
					# add traversal to make sure ID is unique
					allele_id = define_id(ref_string, alt_string, fields[0], ref_pos, traversal_to_string(allele))
					id_to_alleles[allele_id] = (ref_string, alt_string, ref_pos, traversal_to_string(allele), traversal_to_string(ref_allele))
					seen_variants[(ref_string, alt_string, ref_pos)] = allele_id
				id_to_index[allele_id].append(i+1)
				allele_to_ids[i+1].append(allele_id)


		# generate one line per single ID with adjusted genotypes
		for allele_id in id_to_alleles:
			updated_fields = [
				fields[0], # CHROM
				str(id_to_alleles[allele_id][2]), # POS
				fields[2],
				id_to_alleles[allele_id][0], # REF
				id_to_alleles[allele_id][1], # ALT
				'.',
				'PASS',
				'AT=' + id_to_alleles[allele_id][4] + ',' + id_to_alleles[allele_id][3] + ';ID=' + allele_id,
				'GT'
			]
			# update the genotypes
			for genotype in fields[9:]:
				if genotype in ['.', './.']:
					updated_fields.append(genotype)
				else:
					haps = [int(a) if a != '.' else a for a in genotype.split('|')]
					new_genotype = []
					for h in haps:
						if h == '.':
							new_genotype.append('.')
						elif h in id_to_index[allele_id]:
							new_genotype.append('1')
						else:
							new_genotype.append('0')
					updated_fields.append('|'.join(new_genotype))
			biallelic_records.append('\t'.join(updated_fields))
		# generate multiallelic record that is annotated with IDs
		# add ID field to INFO
		info_fields['ID'] = ','.join([':'.join(allele_to_ids[i]) for i in range(1,nr_alleles)])
		fields[7] = info_to_string(info_fields)
		return '\t'.join(fields), biallelic_records


def preprocess_vcf(filename):
	"""
	Read VCF once and store all segment IDs.
	"""
	segments = {}
	with open(filename, 'r') as vcffile:
		for line in vcffile:
			if line.startswith('#'):
				continue
			fields = line.split()
			info_fields = parse_info(fields)
			assert 'LV' in info_fields
#			if int(info_fields['LV']) > 0:
#				# only consider LV==0 variants
#				continue
			assert 'AT' in info_fields
			for node_id in re.split('<|>|,', info_fields['AT']):
				segments[node_id] = True
	return segments
		


if __name__== "__main__":
	parser = argparse.ArgumentParser(prog='annotate_graph.py', description=__doc__)
	parser.add_argument('-vcf', metavar='VCF', required=True, help='VCF file containing top level variants (vcfbub output).')
	parser.add_argument('-gfa', metavar='GFA', required=True, help='graph in GFA format.')
	parser.add_argument('-o', metavar='OUTPREFIX', required=True, help='Prefix of the output files.')
	args = parser.parse_args()

	# stats
	records_read = 0
	multi_written = 0
	bi_written = 0
	
	# read VCF once only to collect the segments needed.
	# this is done to safe space since we do not have to store all graph nodes
	print('Preprocessing...')
	segments = preprocess_vcf(args.vcf)

	# parse the GFA
	print('Reading sequence information from GFA file...')
	gfa = parse_gfa(args.gfa)
	print('Done reading GFA.')

	print('Annotate the VCF file...')
	with open(args.o + '.vcf', 'w') as out_multi, open(args.o + '_biallelic.vcf', 'w') as out_bi:
		for line in open(args.vcf, 'r'):
			if line.startswith('##'):
				out_multi.write(line)
				out_bi.write(line)
				continue
			if line.startswith('#'):
				header_id = '##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">\n'
				out_multi.write(header_id)
				out_bi.write(header_id)
				out_multi.write(line)
				out_bi.write(line)
				continue
			multi_line, bi_lines = decompose(line, gfa)
			if multi_line is None:
				continue
			multi_written += 1
			out_multi.write(multi_line + '\n')
			for b_line in bi_lines:
				bi_written += 1
				out_bi.write(b_line + '\n')
			records_read += 1
			if  records_read % 1000000 == 0:
				print('Processed ' + str(records_read) + ' VCF records.')
	print('Wrote ' + str(multi_written) + ' multi-allelic records.')
	print('Wrote ' + str(bi_written) + ' bi-allelic records.')


