#!/usr/bin/python

import sys
import os
import pandas as pd

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Individual coverage directory: ', sys.argv[1])
print('Combined coverage directory: ', sys.argv[2])
print('Sample: ', sys.argv[3])
print('\n---\n')

# Run combining
with open('full_sample_list.txt') as f:
    sample_list = f.read().splitlines()

df = pd.read_csv(os.path.join(sys.argv[1], sample_list[0] + '.cov.txt'), sep='\t')
df.columns = ['ContigID', 'Avg_fold', 'Length', 'Ref_GC', 'Covered_percent', 'Covered_bases', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev']
contig_list = df.ContigID.values.tolist()
length_list = df.Length.values.tolist()
df = None

cov_dict = {}
cov_var_dict = {}

for sam in sample_list:
	sam_file = os.path.join(sys.argv[1], sam + '.cov.txt')
	df = pd.read_csv(sam_file, sep='\t')
	cov_dict[sam] = df.Avg_fold.values.tolist()
	df['Variance'] = df['Std_Dev'] ** 2
	cov_var_dict[sam + '-var'] = df.Variance.values.tolist()

cov_df = pd.DataFrame(cov_dict)
cov_df['totalAvgDepth'] = cov_df.sum(axis = 1)
cov_var_df = pd.DataFrame(cov_var_dict)

comb_df = pd.concat([cov_df, cov_var_df], axis=1)
comb_df['contigName'] = contig_list
comb_df['contigLen'] = length_list

colnames = list(comb_df.columns)
colnames.remove('contigName')
colnames.remove('totalAvgDepth')
colnames.remove('contigLen')
colnames = list(set(colnames))
new_colnames =  ['contigName', 'contigLen', 'totalAvgDepth'] + colnames
comb_df = comb_df[new_colnames]

comb_df.to_csv(os.path.join(sys.argv[2], sys.argv[3] + '.comb_cov.txt'), sep='\t', header=True, index=False)
