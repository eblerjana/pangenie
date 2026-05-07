import sys

prev_end = 0
prev_chrom = ""

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	fields = line.strip().split()
	chrom = fields[0]
	start = int(fields[1])
	end = start + len(fields[3])
	if (chrom == prev_chrom) and (start < prev_end):
		# record overlaps previous, skip
		continue
	prev_end = end
	prev_chrom = chrom
	print(line.strip())
