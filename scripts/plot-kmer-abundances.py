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

def geometric(p, value):
	return ((1.0-p)**(float(value)) )* p

parser = argparse.ArgumentParser(prog='plot-kmer-abundances.py', description=__doc__)
parser.add_argument('histo', metavar='HISTO', help='.histo file output of PGGTyper.')
parser.add_argument('--max-value', default='10000', metavar='MAX_VALUE', help='max value to plot on x-axis (default: 10000).')
parser.add_argument('--normalize-poisson', default=False, action='store_true', help='plot Poission distributions such that P(k|cn_0) + P(k|cn_1) + P(k|cn_2) = 1')
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
	val_cn0 = geometric(mean_cn0, x_value)
	val_cn1 = poisson(mean_cn1, x_value)
	val_cn2 = poisson(mean_cn2, x_value)
	sum = val_cn0 + val_cn1 + val_cn2 if args.normalize_poisson else 1.0
	if sum == 0.0:
		sum = 1.0
	p_cn0.append(val_cn0 / sum);
	p_cn1.append(val_cn1 / sum);
	p_cn2.append(val_cn2 / sum);

# plot distributions
suffix = '-normalized' if args.normalize_poisson else ''
plt.plot(x_values, normalized_kmer_abundances, '-', label='distribution of kmer counts')
plt.plot(x_values, p_cn0, '-', label='geometric p = '+str(mean_cn0))
plt.plot(x_values, p_cn1, '-', label='poisson mean = '+str(mean_cn1))
plt.plot(x_values, p_cn2, '-', label='poisson mean = '+str(mean_cn2))
plt.axis([0, int(args.max_value), 0,1])
plt.legend(loc='upper right')
plt.xlabel('kmer count')
plt.ylabel('frequency')
plt.savefig(args.histo + suffix + '.pdf')
plt.close()
