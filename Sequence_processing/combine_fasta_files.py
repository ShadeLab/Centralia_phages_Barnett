#!/usr/bin/python

import sys
import os
import pandas as pd
import re
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Sample list: ', sys.argv[1])
print('Annotation directory: ', sys.argv[2])
print('Reference sequences: ', sys.argv[3])
print('\n---\n')

sample_list_file = sys.argv[1]
annotation_dir =  sys.argv[2]
ref_file = sys.argv[3]

with open(sample_list_file, 'r') as f:
	sample_list = [line.rstrip() for line in f]

with open(os.path.join(annotation_dir, 'viral_blast', 'combined_viral_contigs_RefSeq.fasta'), 'w') as combfasta:
	with open(os.path.join(annotation_dir, 'viral_blast', 'combined_viral_contigs.fasta'), 'w') as contigfasta:
		with open(os.path.join(annotation_dir, 'viral_blast', 'contig_name_map.txt'), 'w') as contig_names:
			contig_names.write('SequenceID\tcontig_id\tcontig_number\n')
			contig_num = 0
			for sam in sample_list:
				contig_fasta = SeqIO.parse(os.path.join(annotation_dir, 'virsorter_round2', sam + '.vs2_round2.out', 'final-viral-combined.fa'), 'fasta')
				for record in contig_fasta:
					contig_num += 1
					combfasta.write('>' + 'viral_contig_' + str(contig_num) + '\n')
					combfasta.write(str(record.seq) + '\n')
					contigfasta.write('>' + 'viral_contig_' + str(contig_num) + '\n')
					contigfasta.write(str(record.seq) + '\n')
					contig_names.write(sam + '\t' + record.id + '\t' + str(contig_num) + '\n')
				contig_fasta = None
		print('Number of potentially viral contigs: ' + contig_num)
	RefSeq_fasta = SeqIO.parse(ref_file, 'fasta')
	for record in RefSeq_fasta:
		combfasta.write('>' + record.id + '\n')
		combfasta.write(str(record.seq) + '\n')
	RefSeq_fasta = None
