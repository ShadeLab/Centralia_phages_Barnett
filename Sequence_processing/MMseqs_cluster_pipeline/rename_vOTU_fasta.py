#!/usr/bin/python

import sys
import os
import pandas as pd
import re
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Original vOTU fasta: ', sys.argv[1])
print('New vOTU fasta: ', sys.argv[2])
print('Original to new name map: ', sys.argv[3])
print('\n---\n')

input_fasta = sys.argv[1]
output_fasta =  sys.argv[2]
map_file = sys.argv[3]

with open(map_file, 'w') as outmap:
	outmap.write('vOTU\tShortID\n')
	with open(output_fasta, 'w') as outfasta:
		infasta = SeqIO.parse(input_fasta, 'fasta')
		vOTU_num = 1
		for record in infasta:
			outfasta.write('>' + 'vOTU_' + str(vOTU_num) + '\n')
			outfasta.write(str(record.seq) + '\n')
			outmap.write(record.id + '\t' + 'vOTU_' + str(vOTU_num) + '\n')
			vOTU_num += 1

print('Number of vOTU: ' + str(vOTU_num - 1))
