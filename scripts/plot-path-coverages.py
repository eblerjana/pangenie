import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict



def plot(path_to_fraction, cur_chrom, pdf):
	labels = [p for p in sorted(path_to_fraction.keys())]
	data = [path_to_fraction[p] for p in sorted(path_to_fraction.keys())]

	fig, ax = plt.subplots(figsize=(15,5))
	ax.violinplot(data)
	ax.set_title('Path coverages for chromosome ' + cur_chrom)
	ax.set_xticks([i for i in range(1, len(data)+1)])
	ax.set_xticklabels(labels, rotation=90)
	plt.tight_layout()
	pdf.savefig()
	plt.close()


outname = sys.argv[1]

## Statistics for chr1
#Ranked haplotypes for bubble position: 15951
#path12  0.354839        31

path_to_fraction = defaultdict(list)
cur_chrom = None


with PdfPages(outname) as pdf:

	for line in sys.stdin:
		if line.startswith('#'):
			if cur_chrom:
				# plot violinplot for previous chromosome
				plot(path_to_fraction, cur_chrom, pdf)
			cur_chrom = line.strip().split()[-1]
			path_to_fraction = defaultdict(list)
			continue
		if line.startswith('Ranked'):
			continue
		fields = line.strip().split()
		path = fields[0]
		fraction = float(fields[1])
		path_to_fraction[path].append(fraction)
	plot(path_to_fraction, cur_chrom, pdf)