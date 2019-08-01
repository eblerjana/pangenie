import argparse
import math
import numpy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='plot-kmer-abundances.py', description=__doc__)
parser.add_argument('histo', metavar='HISTO', help='.histo file output of PGGTyper.')
parser.add_argument('--max-value', default='10000', metavar='MAX_VALUE', help='max value to plot on x-axis (default: 10000).')
args = parser.parse_args()

# read the histogram from input file
x_values = []
kmer_abundances = []
for line in open(args.histo, 'r'):
	parsed_line = line.split();
	if parsed_line[0] == "parameters":
		continue
	x_value = int(parsed_line[0])
	y_value = int(parsed_line[1])
	x_values.append(x_value);
	kmer_abundances.append(y_value);
# plot distributions
plt.plot(x_values[1:], kmer_abundances[1:], '-')
plt.xlim(1,int(args.max_value))
plt.ylim(1, max(kmer_abundances)*10)
plt.yscale('log')
plt.title('histogram of kmer abundances')
plt.xlabel('kmer abundance')
plt.ylabel('count')
plt.savefig(args.histo + '.pdf')
plt.close()
