#!/usr/bin/env python
import sys
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import upsetplot
from scipy.stats import pearsonr
from pylab import *
from collections import defaultdict

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import AdaBoostClassifier
from sklearn.svm import SVR
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer



if __name__ == "__main__":
	filename = sys.argv[1]
	outname = sys.argv[2]
#	med_svs = sys.stgv[3]

	df = pd.read_csv(filename, sep='\t')
	df = df.assign(pangenie_mendelian_consistency=lambda d: d['pangenie_mendelian_consistent_trios'] / d['pangenie_considered_trios'])

	# consider allele frequency computed across all genotyped samples
	df = df.assign(ac0_fail = lambda df: df['pangenie-all_allele_freq'] == 0)
	df = df.assign(mendel_fail = lambda df: (df.pangenie_mendelian_consistency < 0.8) & (df['pangenie_considered_trios']>=5))
	df = df.assign(gq_fail = lambda df: df['pangenie-all_GQ>=200'] < 50)
	df = df.assign(self_fail = lambda df: (df['pangenie_self-genotyping_correct [%]'] < 90.0) )
	df = df.assign(nonref_fail = lambda df: ( (df['pangenie_self-genotyping_0/1_typed_0/1'] + df['pangenie_self-genotyping_1/1_typed_1/1'])==0) & ( (df['pangenie_self-genotyping_0/1_typed_0/0'] + df['pangenie_self-genotyping_0/1_typed_1/1'] + df['pangenie_self-genotyping_1/1_typed_0/0'] + df['pangenie_self-genotyping_1/1_typed_0/1']) != 0) )
	df = df.assign(all_pass = lambda df: ~(df.ac0_fail | df.mendel_fail | df.gq_fail | df.nonref_fail | df.self_fail) )
	df = df.assign(negative = lambda df: ~df.ac0_fail & (((1*df.gq_fail) + (1*df.mendel_fail) + (1*df.nonref_fail) + (1*df.self_fail)) >= 3))

	snps = set(id for id in df.variant_id if 'SNV' in id)
	indels = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])<50))
	small_indels = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])<=19))
	midsize_indels = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])>19) and (int(id.split('-')[-1])<50))
	svs = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])>=50))

	small_insertions = set(id for id in small_indels if id.split('-')[2]=="INS")
	small_deletions = set(id for id in small_indels if id.split('-')[2]=="DEL")
	small_complex = set(id for id in small_indels if id.split('-')[2]=="COMPLEX")
	assert len(small_insertions) + len(small_deletions) + len(small_complex) == len(small_indels)

	midsize_insertions = set(id for id in midsize_indels if id.split('-')[2]=="INS")
	midsize_deletions = set(id for id in midsize_indels if id.split('-')[2]=="DEL")
	midsize_complex = set(id for id in midsize_indels if id.split('-')[2]=="COMPLEX")
	assert len(midsize_insertions) + len(midsize_deletions) + len(midsize_complex) == len(midsize_indels)

	large_insertions = set(id for id in svs if id.split('-')[2]=="INS")
	large_deletions = set(id for id in svs if id.split('-')[2]=="DEL")
	large_complex = set(id for id in svs if id.split('-')[2]=="COMPLEX")
	assert len(large_insertions) + len(large_deletions) + len(large_complex) == len(svs)

#	med = set()
#	for line in open(med_svs, 'r'):
#		if not 'allele-' in line:
#			med.add(line.strip())

	# for pangenie, consider allele frequencies computed across unrelated samples (+ panel samples) only, and skip children (unless they are panel samples)
	metrics = ['panel_allele_freq', 'pangenie-unrelated_allele_freq', 'pangenie-unrelated_heterozygosity']

	# features used for regression
	features = [
		'pangenie_self-genotyping_correct [%]',
		'pangenie_self-genotyping_wrong [%]',
		'pangenie_self-genotyping_not_typed [%]',
		'pangenie_self-genotyping_correct',
		'pangenie_self-genotyping_wrong',
		'pangenie_self-genotyping_not_typed',
		'pangenie_self-genotyping_absent_in_truth',
		'pangenie_self-genotyping_0/0_typed_0/0',
		'pangenie_self-genotyping_0/0_typed_0/1',
		'pangenie_self-genotyping_0/0_typed_1/1',
		'pangenie_self-genotyping_0/0_typed_./.',
		'pangenie_self-genotyping_0/1_typed_0/0',
		'pangenie_self-genotyping_0/1_typed_0/1',
		'pangenie_self-genotyping_0/1_typed_1/1',
		'pangenie_self-genotyping_0/1_typed_./.',
		'pangenie_self-genotyping_1/1_typed_0/0',
		'pangenie_self-genotyping_1/1_typed_0/1',
		'pangenie_self-genotyping_1/1_typed_1/1',
		'pangenie_self-genotyping_1/1_typed_./.',
		'panel_allele_freq',
		'panel_alternative_alleles',
		'panel_total_alleles',
		'pangenie-all_alternative_alleles',
		'pangenie-all_total_alleles',
		'pangenie-all_heterozygosity',
		'pangenie-all_heterozygous_genotypes',
		'pangenie-all_total_genotypes',
		'pangenie-all_unique_kmers',
		'pangenie-all_GQ>=200',
		'pangenie-all_allele_freq',
		'pangenie_mendelian_consistent_trios',
		'pangenie_alternative_transmitted',
		'pangenie_considered_trios',
	]

	variant_types = {
		'snps': [snps, ['unfiltered', 'strict'], ['all-regions']],
#		'small_insertions': [small_insertions, ['unfiltered', 'strict'], ['all-regions']],
#		'small_deletions': [small_deletions, ['unfiltered', 'strict'], ['all-regions']],
#		'small_complex': [small_complex, ['unfiltered', 'strict'], ['all-regions']],
#		'midsize_deletions': [midsize_deletions, ['unfiltered', 'strict'], ['all-regions']],
#		'midsize_insertions': [midsize_insertions, ['unfiltered', 'strict'], ['all-regions']],
#		'midsize_complex': [midsize_complex, ['unfiltered', 'strict'], ['all-regions']],
		'indels': [indels, ['unfiltered', 'strict'], ['all-regions']],
		'large_deletions': [large_deletions, ['unfiltered', 'lenient_-0.5', 'strict'], ['all-regions']],
		'large_insertions': [large_insertions, ['unfiltered', 'lenient_-0.5', 'strict'], ['all-regions']],
		'large_complex': [large_complex, ['unfiltered', 'lenient_-0.5', 'strict'], ['all-regions']]
#		'giab_med_svs': [med, ['unfiltered', 'lenient_-0.5', 'strict']]
	}

	df['score_SVR'] = np.nan

	for name, info in variant_types.items():
		if 'large' in name:
			df_sub = df[df.variant_id.isin(info[0]) & ~df.ac0_fail].copy()
			df_sub.sort_values(by=['variant_id'], inplace=True)
			df_sub.set_index('variant_id')

			# impute missing values
			imp = IterativeImputer(max_iter=10, verbose=0, random_state=0)
			imp.fit(df_sub[features])
			imputed_df = pd.DataFrame(imp.transform(df_sub[features]), columns=df_sub[features].columns, index=df_sub.index)
			df_sub.update(imputed_df, overwrite=True)

			# Fit transform using all data points (labeled + unlabeled)
			df_sub = df_sub.assign(target = lambda df: (-1*df.negative) + (1*df.all_pass))
			x = df_sub.loc[:,features].values
			scaler = StandardScaler()
			scaler.fit_transform(x)

			# Train model only on labeled points
			autosomal_ids = [i for i in info[0] if not 'chrX' in i]
			df_labeled = df_sub[df.variant_id.isin(autosomal_ids) & (df_sub.target != 0)]
			x = scaler.transform(df_labeled.loc[:,features].values)
			y = df_labeled.loc[:,['target']].values
			print('Training regression model')
			regressor = SVR(C=50)
#			regressor = GradientBoostingRegressor(n_estimators=20, random_state=0)
			regressor.fit(x,y.ravel())

			# Apply to unlabeled data
			y_pred = regressor.predict(scaler.transform(df_sub.loc[:,features].values))
			column_label = "score_" + name
			df_score = pd.DataFrame({"variant_id":df_sub.variant_id, column_label:y_pred})

			# Add column with variant specific scores to table
			df = df.merge(df_score, on='variant_id', how='left')

		for filter in info[1]:
			variant = info[0]
			for region in info[2]:
				# collect all variants in the region
				if region == "all-regions":
					ids_region = variant
				elif region == "repeat-regions":
#					ids_region = set([i for i in info[0] if (df.loc[i, 'ucsc_overlaps'] > 0)])
					ids_region = set(id for id in df[df["ucsc-repeats_overlaps"]>=0.5].variant_id if id in variant)
				elif region == "nonrepeat-regions":
#					ids_region = set([i for i in info[0] if (df.loc[i, 'ucsc_overlaps'] <= 0)])
					ids_region = set(id for id in df[df["ucsc-repeats_overlaps"]<0.5].variant_id if id in variant)

				# create plots
				with PdfPages(outname + '_' + name + '_' + filter + '_' + region + '.pdf') as pdf:
					ids_autosomes = set([i for i in ids_region if not 'chrX' in i])
				
					if filter == 'unfiltered':
						df_sub = df[df.variant_id.isin(ids_region)]
						df_sub_autosomes = df[df.variant_id.isin(ids_autosomes)]
					elif filter.startswith('lenient'):
						score_column = 'score_' + name
						plt.figure()
						fig, ax = plt.subplots()
						df[df.variant_id.isin(variant) & df.negative].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='red')
						df[df.variant_id.isin(variant) & df.all_pass].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='blue')
						df[df.variant_id.isin(variant) & ~df.negative & ~df.all_pass & ~df.ac0_fail].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='grey')
						ax.set_yscale('log')
						pdf.savefig()
						plt.close()
	
						score_cutoff = float(filter.split('_')[1])
						df_sub = df[df.variant_id.isin(ids_region) & ( df.all_pass | ((~df.ac0_fail) & (df[score_column] > score_cutoff)) ) ]
						df_sub_autosomes = df[df.variant_id.isin(ids_autosomes) & ( df.all_pass | ((~df.ac0_fail) & (df[score_column] > score_cutoff)) ) ]
					else:
						assert(filter == 'strict')
						df_sub = df[df.variant_id.isin(ids_region) & df.all_pass]
						df_sub_autosomes = df[df.variant_id.isin(ids_autosomes) & df.all_pass]

					print('  variant count ' + filter + ' ' + region + ' ' + name + ':', len(df_sub))

					# create upset plot with filters
					filter_counts = df_sub.groupby(by=['ac0_fail','mendel_fail','gq_fail','nonref_fail', 'self_fail']).size()
					plt.figure()
					upsetplot.plot(filter_counts, sort_by='cardinality')
					pdf.savefig()
					plt.close()


					for metric in ['pangenie_mendelian_consistency', 'pangenie-unrelated_allele_freq']:
						plt.figure()
						fig, ax = plt.subplots()
						df_sub.hist(ax=ax, column=metric, bins=64, bottom = 0.1)
						ax.set_yscale('log')
						fig.suptitle('min={}, max={}, mean={}, median={}'.format(df_sub[metric].min(),df_sub[metric].max(),df_sub[metric].mean(),df_sub[metric].median()) )
						pdf.savefig()
						plt.close()

					# plot panel allele frequencies for variants typed with AF=0
					df_absent = df_sub[df_sub['pangenie-all_allele_freq'] == 0.0]
					print(' typed AF=0 (allowing missing genotypes) ' + filter + ' ' + region + ' ' + name + ':', len(df_absent))
					print(' typed 0/0 in all samples ' + filter + ' ' + region + ' ' + name + ':', len(df_absent[df_absent['pangenie-all_unknown_alleles'] == 0]))
					plt.figure()
					fig, ax = plt.subplots()
					df_absent.hist(ax=ax, column='panel_allele_freq', bins=64, bottom = 0.1)
					ax.set_yscale('log')
					fig.suptitle('panel allele freq. for variants typed with AF=0')
#					fig.suptitle('min={}, max={}, mean={}, median={}'.format(df_absent['panel_allele_freq'].min(),df_absent['panel_allele_freq'].max(),df_absent['panel_allele_freq'].mean(),df_absent['panel_allele_freq'].median()) )
					pdf.savefig()
					plt.close()


					# heatmaps
					for i in range(len(metrics)):
						for j in range(i+1, len(metrics)):
								plt.figure()
								fig, ax = plt.subplots()
								x_values = []
								y_values = []
								for l,f in zip(df_sub_autosomes[metrics[i]], df_sub_autosomes[metrics[j]]):
									if not pd.isnull(f) and not pd.isnull(l):
										x_values.append(l)
										y_values.append(f)
								assert len(x_values) == len(y_values)
								if len(x_values) == 0:
									continue
	#							cmap = cm.get_cmap('Greys', 6) 
	#							joint_kws=dict(gridsize=35, cmap=cmap)
								joint_kws=dict(gridsize=35, cmap="hot_r")
								ax = sns.jointplot(x=x_values, y=y_values, xlim=(-0.05, 1.05), ylim=(-0.05,1.05), bins='log', kind='hex', joint_kws=joint_kws, marginal_ticks=True, color="red")
								ax.set_axis_labels(metrics[i], metrics[j])
								if 'allele_freq' in metrics[i] and 'allele_freq' in metrics[j]:
									pearson_corr, p_value = pearsonr(x_values, y_values)
									ax.fig.suptitle(name + " (r=" + str(pearson_corr) + ")")
									print('  pearson correlation ' + name + ' ' + filter + ' ' + region + ':', metrics[i], metrics[j], pearson_corr)
								if 'pangenie-unrelated_allele_freq' in [metrics[i], metrics[j]] and 'pangenie-unrelated_heterozygosity' in [metrics[i], metrics[j]]:
									# plot theoretical line
									t = np.arange(0.0, 1.01, 0.01)
									s = [2*j*(1.0-j) for j in t]
									ax.fig.suptitle(name)
									ax.ax_joint.plot(t,s,'r-')
								plt.colorbar()
								plt.tight_layout()
								pdf.savefig()
								plt.close()
			
					# boxplots
					pop_metrics = ['pangenie-unrelated_allele_freq', 'panel_allele_freq']
					length_intervals = []
					n_intervals = 10
					previous = -0.1
					length_intervals = []
					for i in range(1, n_intervals + 1):
						length_intervals.append([previous, i/n_intervals + 0.0001])
						previous = i/n_intervals
					print(length_intervals)
					length_afs = [[] for i in range(n_intervals)]
					for l,f in zip(df_sub_autosomes[pop_metrics[0]], df_sub_autosomes[pop_metrics[1]]):
						if not pd.isnull(l):
							for i,interval in enumerate(length_intervals):
								if interval[0] < f <= interval[1]:
									length_afs[i].append(l)
									break
	
					print(sum([len(b) for b in length_afs]))
					# create boxplots
					fig = plt.figure()
					ax = plt.axes()		
					bp = ax.boxplot(length_afs)
					length_labels = []
					for i,x in enumerate(length_intervals):
						label = '<=' + str(i+1) + '/' + str(n_intervals)
						length_labels.append(label)
	#				ax.set_yticks([0] + [i[1] for i in length_intervals])
	#				ax.set_yticklabels(["0"] + [str(i) + '/' + str(n_intervals) for i in range(1,n_intervals+1)])
					ax.set_xticklabels(length_labels, rotation=45)
			#		ax.set_yscale("log")
					plt.tight_layout()
					pdf.savefig()
					plt.close()


					# plot distribution of AF=0 variants along the chromosomes
					# do this for large variants only
					if (filter == 'unfiltered') and ('large' in name):
						chromosome_to_total = defaultdict(list)
						chromosome_to_absent = defaultdict(list)

						for var_id in df_sub.variant_id:
							fields = var_id.split('-')
							chromosome_to_total[fields[0]].append(int(fields[1]))
						for var_id in df_absent.variant_id:
							fields = var_id.split('-')
							chromosome_to_absent[fields[0]].append(int(fields[1]))

						# create histogram for each chromosome
						for chrom in sorted(chromosome_to_total.keys()):
							fig = plt.figure(figsize=(80,10))
							ax = plt.axes()
							fig.suptitle(chrom + ' ' + name)

							_, bins, _ = plt.hist(chromosome_to_total[chrom], bins=500, alpha=0.5, label="total alleles")
							_ = plt.hist(chromosome_to_absent[chrom], bins=bins, alpha=0.5, label="AF=0 alleles")

							plt.xlabel('genomic position (bp)')
							plt.ylabel('allele count')
#							plt.yscale('log', nonposy='clip')
							plt.legend(loc='upper right')
							plt.tight_layout()
							pdf.savefig()
							plt.close()






	for variant in variant_types.keys():
		if 'large' in variant:
			column_name = 'score_' + variant
			df['score_SVR'].fillna(df[column_name], inplace=True)
	
	df['confidence_level'] = 0
	df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<-0.5)), 1, inplace=True )
	df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<0.0)), 2, inplace=True )
	df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<0.5)), 3, inplace=True )
	df.confidence_level.where(~df.all_pass, 4, inplace=True)

	df = df.assign(is_snp = lambda df: df.variant_id.isin(snps))
	df = df.assign(is_small_insertion = lambda df: df.variant_id.isin(small_insertions))
	df = df.assign(is_small_deletion = lambda df: df.variant_id.isin(small_deletions))
	df = df.assign(is_small_complex = lambda df: df.variant_id.isin(small_complex))
	df = df.assign(is_midsize_insertion = lambda df: df.variant_id.isin(midsize_insertions))
	df = df.assign(is_midsize_deletion = lambda df: df.variant_id.isin(midsize_deletions))
	df = df.assign(is_midsize_complex = lambda df: df.variant_id.isin(midsize_complex))
	df = df.assign(is_large_insertion = lambda df: df.variant_id.isin(large_insertions))
	df = df.assign(is_large_deletion = lambda df: df.variant_id.isin(large_deletions))
	df = df.assign(is_large_complex = lambda df: df.variant_id.isin(large_complex))

	header = ["variant_id", "ac0_fail", "mendel_fail", "self_fail", "gq_fail", "nonref_fail", "score_SVR", "all_pass", "confidence_level", "is_snp", "is_small_insertion", "is_small_deletion", "is_small_complex", "is_midsize_insertion", "is_midsize_deletion", "is_midsize_complex", "is_large_insertion", "is_large_deletion", "is_large_complex"]
	df.to_csv(outname + '_filters.tsv', columns=header, sep='\t', index=False, na_rep='nan')

