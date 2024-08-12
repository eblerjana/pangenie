import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict



def plot(path_to_fraction, cur_chrom, pdf, caption):
	labels = [p for p in sorted(path_to_fraction.keys())]
	data = [path_to_fraction[p] for p in sorted(path_to_fraction.keys())]

	fig, ax = plt.subplots(figsize=(15,5))
	ax.violinplot(data)
	ax.set_title(caption + cur_chrom)
	ax.set_xticks([i for i in range(1, len(data)+1)])
	ax.set_xticklabels(labels, rotation=90)
	plt.tight_layout()
	pdf.savefig()
	plt.close()



def process_bubble(lines):
	if not lines:
			return {}, False
	path_to_fraction = {}
	allele_to_fraction = {}
	allele_to_count = defaultdict(lambda: 0)
	nr_paths = 0
	consider_bubble = False

	for line in lines:
		nr_paths += 1
		fields = line.strip().split()
		path = fields[0]
		allele = int(fields[2])
		fraction = float(fields[1])
		path_to_fraction[path] = fraction
		allele_to_fraction[allele] = fraction
		allele_to_count[allele] += 1
	
	rare_alleles = []
	for allele in allele_to_count.keys():
		a_count = allele_to_count[allele]
		a_freq = allele_to_fraction[allele]
		af = a_count / nr_paths
		if (af <= 0.05) and (a_freq >= 0.9):
			# rare variant that is well
			# supported by kmers
			rare_alleles.append(allele)
	
	# make sure there is a single high scoring rare allele
	if len(rare_alleles) == 1:
		for allele in allele_to_count.keys():
			if allele == rare_alleles[0]:
				continue
			if allele_to_fraction[allele] > 0.3:
				# found other high scoring allele
				break
		consider_bubble = True

	if consider_bubble:
		print("-------------------------------------------")
		for l in lines:
			print(l)
		print("-------------------------------------------")
	return path_to_fraction, consider_bubble
		




def add_stats(lines, path_to_fraction_all, path_to_fraction_present):
	path_to_fraction, consider_bubble = process_bubble(lines)
	for k,v in path_to_fraction.items():
		path_to_fraction_all[k].append(v)
		if consider_bubble:
			path_to_fraction_present[k].append(v)


if __name__ == "__main__":

	outname = sys.argv[1]
	path_to_fraction_all = defaultdict(list)
	path_to_fraction_present = defaultdict(list)
	cur_chrom = None

	lines = []


	with PdfPages(outname) as pdf:

		for line in sys.stdin:
			if line.startswith('#'):
				add_stats(lines, path_to_fraction_all, path_to_fraction_present)
				if cur_chrom:
					# plot violinplot for previous chromosome
					plot(path_to_fraction_all, cur_chrom, pdf, 'Path coverages (all bubbles) for chromosome ')
					plot(path_to_fraction_present, cur_chrom, pdf, 'Path coverages (non 0/0 bubbles) for chromosome ')
				cur_chrom = line.strip().split()[-1]
				path_to_fraction_present = defaultdict(list)
				path_to_fraction_all = defaultdict(list)
				lines = []
				continue
			if line.startswith('Ranked'):
				# process previously read bubble
				add_stats(lines, path_to_fraction_all, path_to_fraction_present)
				lines = []
				continue
			lines.append(line)


		# process the last record
		add_stats(lines, path_to_fraction_all, path_to_fraction_present)

		plot(path_to_fraction_all, cur_chrom, pdf, 'Path coverages (all bubbles) for chromosome ')
		plot(path_to_fraction_present, cur_chrom, pdf, 'Path coverages (non 0/0 bubbles) for chromosome ')