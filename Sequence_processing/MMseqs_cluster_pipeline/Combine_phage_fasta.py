#!/usr/bin/python

import sys
import os
import pandas as pd
import re
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Viral contig list: ', sys.argv[1])
print('Phage directory: ', sys.argv[2])
print('Contig directory: ', sys.argv[3])
print('\n---\n')

viral_contig_list_file = sys.argv[1]
phage_contig_dir =  sys.argv[2]
contig_dir = sys.argv[3]

viral_contig_list = pd.read_csv(viral_contig_list_file, sep='\t')
viral_contig_list = viral_contig_list[viral_contig_list['is_phage'] == 1]

sample_list = list(set(viral_contig_list['SequenceID']))

with open(os.path.join(phage_contig_dir, 'FullSet_viral_contigs.fasta'), 'w') as out_fasta:
	for sam in sample_list:
		viral_contigs = viral_contig_list[viral_contig_list['SequenceID'] == sam]
		print('In ' + sam + ' there should be ' + str(len(viral_contigs['viralID'])) + ' viral contigs')
		contig_fasta = SeqIO.parse(os.path.join(contig_dir, sam + '.vs2_round2.out', 'final-viral-combined.fa'), 'fasta')
		viral_count = 0
		for record in contig_fasta:
			viral_id = record.id.replace('|', '_')
			if viral_id in list(viral_contigs['viralID']):
				viral_count += 1
				out_fasta.write('>' + sam + '__' + viral_id + '\n')
				out_fasta.write(str(record.seq) + '\n')
		contig_fasta = None
		print('We actually found ' + str(viral_count) + ' viral contigs in the fasta file\n')
	viral_contigs = None
