#!/usr/bin/python

import sys
import os
import pandas as pd

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('CRISPR file: ', sys.argv[1])
print('Output CRISPR fasta: ', sys.argv[2])
print('\n---\n')

crispr_file = sys.argv[1]
crispr_fasta = sys.argv[2]

# Make fasta file from CRISPR dataframe
df = pd.read_csv(crispr_file, sep='\t')

with open(crispr_fasta, 'w') as outfile:
	for index, row in df.iterrows():
		outfile.write('>' + '__'.join([row['SequenceID'], row['ContigID'], 'Array_' + str(row['CRISPR_ID']), 'Spacer_' + str(row['Spacer_ID'])]) + '\n')
		outfile.write(row['Spacer_seq'] + '\n')
