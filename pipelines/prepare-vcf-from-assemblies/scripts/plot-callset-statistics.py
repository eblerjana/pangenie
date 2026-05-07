import sys
import argparse
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt 
from collections import namedtuple

SampleStats = namedtuple('SampleStats', 'snps insertions deletions other het_hom ti_tv')

class Statistics:
	def __init__(self, vcfstats_file):
		self._statistics = {}
		current_sample = None
		nr_snps = 0
		nr_insertions = 0
		nr_deletions = 0
		nr_other = 0
		het_hom = 0.0
		ti_tv = 0.0
		for line in open(vcfstats_file, 'r'):
			fields = line.split()
			if line.startswith('Sample Name:'):
				if current_sample is not None:
					# store stats of previous sample
					self._statistics[current_sample] = SampleStats(nr_snps, nr_insertions, nr_deletions, nr_other, het_hom, ti_tv)
					current_sample = 0
					nr_snps = 0
					nr_insertions = 0
					nr_deletions = 0
					nr_other = 0
					het_hom = 0.0
					ti_tv = 0.0
				current_sample = fields[-1]
			if line.startswith('SNPs'):
				nr_snps = int(fields[-1])
			if line.startswith('MNPs'):
				nr_other += int(fields[-1])
			if line.startswith('Insertions'):
				nr_insertions = int(fields[-1])
			if line.startswith('Deletions'):
				nr_deletions = int(fields[-1])
			if line.startswith('Indels'):
				nr_other += int(fields[-1])
			if line.startswith('Total Het/Hom ratio'):
				het_hom = float(fields[-2])
			if line.startswith('SNP Transitions/Transversions'):
				ti_tv = float(fields[-2])
		# add last sample
		self._statistics[current_sample] = SampleStats(nr_snps, nr_insertions, nr_deletions, nr_other, het_hom, ti_tv)

	def nr_snps(self, sample):
		return self._statistics[sample].snps

	def nr_insertions(self, sample):
		return self._statistics[sample].insertions

	def nr_deletions(self, sample):
		return self._statistics[sample].deletions

	def nr_others(self, sample):
		return self._statistics[sample].other

	def het_hom(self, sample):
		return self._statistics[sample].het_hom

	def ti_tv(self, sample):
		return self._statistics[sample].ti_tv

	def samples(self):
		return [s for s in self._statistics.keys()]


def run_length(minlen, outname):
	x_values = [i for i in range(-minlen, minlen+1)]
	y_values = [0] * len(x_values)

	for line in sys.stdin:
		if line.startswith('LENGTH'):
			continue
		# histogram line
		fields = line.split()
		assert len(fields) == 3
		length = int(fields[0])
		if length == 0:
			continue
		if abs(length) > minlen:
			continue
		y_values[length + minlen] = int(fields[1])

	# plotting
#	plt.figure(figsize=(20, 15))
	plt.plot(x_values, y_values, linestyle='-')
	plt.xlabel('indel length')
	plt.ylabel('Count')
	plt.xlim(-minlen, minlen)
	plt.yscale('symlog')
	plt.xscale('symlog')
	plt.savefig(outname)

def run_ratios(input, outname, samples):
	# read from rtg vcfstats output
	statistics = Statistics(input)
	if samples == []:
		samples = statistics.samples()
	print(samples, 'samples')
	values = [(statistics.het_hom(s), s) for s in samples]
	# sort by het/hom ratio
	values.sort(key=lambda tup: tup[0])
	tick_values = [i for i in range(len(samples))]
	tick_labels = [s for (i,s) in values]
	x_values = [i for (i,s) in values]

	plt.hlines(y=tick_values, xmin=0, xmax=x_values, color='#570502')
	plt.plot(x_values, tick_values, 'o', color='#570502')
	plt.yticks(tick_values,tick_labels, rotation=45)
	plt.title("het/hom ratios")
	plt.savefig(outname)

def run_numbers(input, outname, samples):
	statistics = Statistics(input)
	if samples == []:
		samples = statistics.samples()
	x_values = [i for i in range(len(samples))]
	x_labels = [s for s in samples]
	snps = [statistics.nr_snps(s) for s in samples]
	insertions = [statistics.nr_insertions(s) for s in samples]
	deletions = [statistics.nr_deletions(s) for s in samples]
	others = [statistics.nr_others(s) for s in samples]

	plt.bar(x_values, snps, label='SNPs', color='#900D04')
	previous = snps
	plt.bar(x_values, insertions, bottom=previous, label='Insertions', color='#F5A66D')
	previous = [i+j for i,j in zip(previous, insertions)]
	plt.bar(x_values, deletions, bottom=previous, label='Deletions', color='#E97303')
	previous = [i+j for i,j in zip(previous, deletions)]
	plt.bar(x_values, others, bottom=previous, label='Others', color='#570502')

	plt.xticks(x_values, x_labels, rotation=45)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc='lower right')
	plt.tight_layout()
	plt.savefig(outname)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot-callset-statistics.py', description=__doc__)
	subparsers = parser.add_subparsers(dest="subparser_name")
	parser_length_hist = subparsers.add_parser('length', help='create indel length histogram.')
	parser_length_hist.add_argument('outname', metavar='OUTNAME', help='name of the output file.')
	parser_length_hist.add_argument('minlen', type=int, metavar='MINLEN', help='minimum indel length.')

	parser_ratios = subparsers.add_parser('ratios', help='plot het/hom ratios of all samples.')
	parser_ratios.add_argument('statistics', metavar='STATISTICS', help='output of rtg vcfstats.')
	parser_ratios.add_argument('outname', metavar='OUTNAME', help='name of the output file.')
	parser_ratios.add_argument('--samples', default=None, help='only consider these samples (given as a comma separated list).')

	parser_numbers = subparsers.add_parser('numbers', help='plot numbers of different variant types.')
	parser_numbers.add_argument('statistics', metavar='STATISTICS', help='output of rtg vcfstats.')
	parser_numbers.add_argument('outname', metavar='OUTNAME', help='name of the output file.')
	parser_numbers.add_argument('--samples', default=None, help='only consider these samples (given as a comma separated list).')
	args = parser.parse_args()

	if args.subparser_name == 'length':
		run_length(args.minlen, args.outname)
	if args.subparser_name == 'ratios':
		samples = args.samples.split(',') if args.samples is not None else []
		print("samples", samples)
		run_ratios(args.statistics, args.outname, samples)
	if args.subparser_name == 'numbers':
		samples = args.samples.split(',') if args.samples is not None else []
		run_numbers(args.statistics, args.outname, samples)
