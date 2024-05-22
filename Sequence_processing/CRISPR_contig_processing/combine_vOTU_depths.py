#!/usr/bin/python

import sys
import os
import pandas as pd

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Individual coverage directory: ', sys.argv[1])
print('Combined coverage directory: ', sys.argv[2])
print('\n---\n')

# Run combining
with open('full_sample_list.txt') as f:
    sample_list = f.read().splitlines()

comb_df = pd.DataFrame()

for sam in sample_list:
	sam_file = os.path.join(sys.argv[1], sam + '.cov.txt')
	df = pd.read_csv(sam_file, sep='\t')
	df.columns = ['ContigID', 'Avg_fold', 'Length', 'Ref_GC', 'Covered_percent', 'Covered_bases', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev']
	df['Variance'] = df['Std_Dev'] ** 2
	df['SequenceID'] = sam
	comb_df = pd.concat([comb_df, df], ignore_index=True)
	df = None

comb_df = comb_df[['SequenceID', 'ContigID', 'Avg_fold', 'Length', 'Ref_GC', 'Covered_percent', 'Covered_bases', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev', 'Variance']]

comb_df.to_csv(os.path.join(sys.argv[2], 'vOTU.comb_cov.txt'), sep='\t', header=True, index=False)
