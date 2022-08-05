import sys
from collections import defaultdict

# list of filters
var_ids_file = sys.argv[1]
# level
level = None
assert sys.argv[2] in ['unfiltered', 'strict', 'lenient']
if sys.argv[2] == "unfiltered":
	level = 0
elif sys.argv[2] == "lenient":
	level = 1
else:
	level = 4

var_ids = defaultdict(lambda: False)

header = None
for line in open(var_ids_file, 'r'):
	if "variant_id" in line:
		header = line.strip().split()
		continue
	fields = line.split()
	index_id = header.index('variant_id')
	index_conf = header.index('confidence_level')
	id = fields[index_id]
	conf = int(fields[index_conf])
	if conf >= level:
		var_ids[id] = True

sys.stderr.write('Read ' + str(len(var_ids)) + ' IDs from input list.\n')

total_lines = 0
skipped_lines = 0
for line in sys.stdin:
	if line.startswith('#'):
		print(line[:-1])
		continue
	fields = line.split()
	info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
	assert 'ID' in info_fields
	assert len(info_fields['ID'].split(',')) == 1
	variant_id = info_fields['ID'].split(',')[0]
	total_lines += 1
	if not var_ids[variant_id]:
		skipped_lines += 1
		continue
	print(line[:-1])
sys.stderr.write('Skipped ' + str(skipped_lines) + ' of ' + str(total_lines) + ' from input vcf.\n')	
