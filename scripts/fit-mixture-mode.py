import argparse
import matplotlib.pyplot as plt
import math
import numpy as np

from pomegranate import *

def poisson(mean, value):
	s = 0.0
	for i in range(1, value+1):
		s += math.log(i)
	log_val = -mean + float(value)*math.log(mean) - s
	return math.exp(log_val)

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
#		if x_value < 6:
#			continue
#		if x_value == 6:
#			for i in range(0,6):
#				all_values.extend([i] * y_value)
#				sum += y_value
		all_values.extend([x_value] * y_value)
		sum += y_value
# fit mixture model
X = np.array([all_values]).T.copy()
print('Fit Mixture Model ...')
pois_0 = PoissonDistribution(mean_cn0)
pois_1 = PoissonDistribution(mean_cn1)
pois_2 = PoissonDistribution(mean_cn2)
model = GeneralMixtureModel([pois_0, pois_1, pois_2])
print(model)
model.fit(X)
print(model)
# mixture weights computed
weight_cn0 = math.exp(model.weights[0])
weight_cn1 = math.exp(model.weights[1])
weight_cn2 = math.exp(model.weights[2])
updated_mean_cn0 = model.distributions[0].parameters[0]
updated_mean_cn1 = model.distributions[1].parameters[0]
updated_mean_cn2 = model.distributions[2].parameters[0]
# print summary
print('########## Summary ##########')
print('Parameters computed by PGGTyper:\tCN0:' + str(mean_cn0) + ' CN1: ' + str(mean_cn1) + ' CN2: ' + str(mean_cn2))
print('Updated parameters of the fitted Poisson mixture:\tCN0: ' + str(updated_mean_cn0) + '\tCN1: ' + str(updated_mean_cn1) + ' CN2: ' + str(updated_mean_cn2))
print('Weights of the fitted Poisson mixture:\tCN0: ' + str(weight_cn0) + '\tCN1: ' + str(weight_cn1) + ' CN2: ' + str(weight_cn2))
# plot the results
p_cn0 = []
p_cn1 = []
p_cn2 = []
for x_value in x_values:
	val_cn0 = weight_cn0 * poisson(updated_mean_cn0, x_value) * sum
	val_cn1 = weight_cn1 * poisson(updated_mean_cn1, x_value) * sum
	val_cn2 = weight_cn2 * poisson(updated_mean_cn2, x_value) * sum
	p_cn0.append(val_cn0)
	p_cn1.append(val_cn1)
	p_cn2.append(val_cn2)
# plot distributions
plt.plot(x_values, kmer_abundances, '-', label='distribution of kmer counts')
plt.plot(x_values, p_cn0, '-', label='Poisson(' + str(updated_mean_cn0) + ')')
plt.plot(x_values, p_cn1, '-', label='Poisson(' + str(updated_mean_cn1) + ')')
plt.plot(x_values, p_cn2, '-', label='Poisson(' + str(updated_mean_cn2) + ')')
plt.xlim(0,int(args.max_value))
plt.legend(loc='upper right')
plt.xlabel('kmer count')
plt.ylabel('number of kmers')
plt.savefig(args.histo + '.pdf')
plt.close()
