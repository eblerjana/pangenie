import argparse
import matplotlib.pyplot as plt
import math
import numpy as np

def poisson(mean, value):
	s = 0.0
	for i in range(1, value+1):
		s += math.log(i)
	log_val = -mean + float(value)*math.log(mean) - s
	return math.exp(log_val)

def geometric(p, value):
	return ((1.0-p)**(float(value)) )* p

parser = argparse.ArgumentParser(prog='fit-mixture-model.py', description=__doc__)
parser.add_argument('histo', metavar='HISTO', help='.histo file output of PGGTyper.')
parser.add_argument('--max-value', default='10000', metavar='MAX_VALUE', help='max value to plot on x-axis (default: 10000).')
args = parser.parse_args()

# read the histogram from input file
print('Read data from file: ' + args.histo + '...')
x_values = []
kmer_abundances = []
all_values = []
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
		all_values.extend([x_value] * y_value)
		sum += y_value

# plot the results
p_cn0 = []
p_cn1 = []
p_cn2 = []
mixture = []
weight_cn0 = 0.7
weight_cn1 = 0.3 / 11.0
weight_cn2 = (0.3 * 10) / 11.0 

for x_value in x_values:
	val_cn0 = weight_cn0 * geometric(mean_cn0, x_value) * sum
	val_cn1 = weight_cn1 * poisson(mean_cn1, x_value) * sum
	val_cn2 = weight_cn2 * poisson(mean_cn2, x_value) * sum
	p_cn0.append(val_cn0)
	p_cn1.append(val_cn1)
	p_cn2.append(val_cn2)
	mixture.append(val_cn0 + val_cn1 + val_cn2)
# plot distributions
plt.semilogy(x_values[1:], kmer_abundances[1:], '-', label='distribution of kmer counts')
plt.semilogy(x_values[1:], p_cn0[1:], '--', label='Geom(' + str(mean_cn0) + ')')
plt.semilogy(x_values[1:], p_cn1[1:], '--', label='Poisson(' + str(mean_cn1) + ')')
plt.semilogy(x_values[1:], p_cn2[1:], '--', label='Poisson(' + str(mean_cn2) + ')')
plt.semilogy(x_values[1:], mixture[1:], '-', label='Mixture')
plt.xlim(1,int(args.max_value))
plt.ylim(1,max(kmer_abundances)*10)
plt.legend(loc='upper right')
plt.xlabel('kmer count')
plt.ylabel('number of kmers')
plt.savefig(args.histo + '.pdf')
plt.close()
