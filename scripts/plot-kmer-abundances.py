import argparse
import math
import numpy
import matplotlib.pyplot as plt

def poisson(mean, value):
	s = 0.0
	for i in range(1, value+1):
		s += math.log(i)
	log_val = -mean + float(value)*math.log(mean) - s
	return math.exp(log_val)

parser = argparse.ArgumentParser(prog='plot-kmer-abundances.py', description=__doc__)
parser.add_argument('histo', metavar='HISTO', help='.histo file output of PGGTyper.')
parser.add_argument('--max-value', default='10000', metavar='MAX_VALUE', help='max value to plot on x-axis (default: 10000).')
args = parser.parse_args()

# read the histogram from input file
x_values = []
kmer_abundances = []
sum = 0.0
mean_cn0 = 0
mean_cn1 = 0
mean_cn2 = 0
for line in open(args.histo, 'r'):
	parsed_line = line.split();
	if parsed_line[0] == 'parameters':
		mean_cn0 = float(parsed_line[1])
		mean_cn1 = float(parsed_line[2])
		mean_cn2 = float(parsed_line[3])
	else:
		x_value = int(parsed_line[0])
		y_value = int(parsed_line[1])
		x_values.append(x_value);
		kmer_abundances.append(y_value);
		sum += y_value
# normalize histogram
normalized_kmer_abundances = [ y/sum for y in kmer_abundances]
# compute poisson probabilities
p_cn0 = []
p_cn1 = []
p_cn2 = []
for x_value in x_values:
	p_cn0.append(poisson(mean_cn0, x_value))
	p_cn1.append(poisson(mean_cn1, x_value))
	p_cn2.append(poisson(mean_cn2, x_value))

# plot distributions
plt.plot(x_values, normalized_kmer_abundances, '-', label='distribution of kmer counts')
plt.plot(x_values, p_cn0, '-', label='poisson mean '+str(mean_cn0))
plt.plot(x_values, p_cn1, '-', label='poisson mean '+str(mean_cn1))
plt.plot(x_values, p_cn2, '-', label='poisson mean '+str(mean_cn2))
plt.axis([0, int(args.max_value), 0,1])
plt.legend(loc='upper right')
plt.xlabel('kmer count')
plt.ylabel('frequency')
plt.savefig(args.histo + '.pdf')
plt.close()
