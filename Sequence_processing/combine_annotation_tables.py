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
print('\n---\n')

sample_list_file = sys.argv[1]
annotation_dir =  sys.argv[2]

with open(sample_list_file, 'r') as f:
	sample_list = [line.rstrip() for line in f]

############################################################
#################### Virsorter2 scores  ####################
############################################################

print('Getting virsorter2 scores')
VS2_scores_df = pd.DataFrame()
for sampleID in sample_list:
	sub_annotation = pd.read_csv(os.path.join(annotation_dir, 'virsorter_round2', sampleID + '.vs2_round2.out', 'final-viral-score.tsv'), sep='\t')
	sub_annotation['SequenceID'] = sampleID
	VS2_scores_df = VS2_scores_df.append(sub_annotation, ignore_index=True)
	sub_annotation = None


print('Number of viral containing contigs ' + str(len(VS2_scores_df.index)))

VS2_scores_df.to_csv(os.path.join(annotation_dir, 'virsorter2_scores.txt'), header=True, index=False, sep='\t')

VS2_scores_df = None

############################################################
################# Contamination scores  ####################
############################################################

print('Getting CheckV contamination scores')
checkV_df = pd.DataFrame()

print('Getting checkV scores')
checkV_df = pd.DataFrame()
for sampleID in sample_list:
        sub_annotation = pd.read_csv(os.path.join(annotation_dir, 'checkV_output', sampleID + '.checkV.out', 'contamination.tsv'), sep='\t')
        sub_annotation['SequenceID'] = sampleID
        checkV_df = checkV_df.append(sub_annotation, ignore_index=True)
        sub_annotation = None


print('Number of CheckV contigs ' + str(len(checkV_df.index)))

checkV_df.to_csv(os.path.join(annotation_dir, 'checkV_contamination.txt'), header=True, index=False, sep='\t')

checkV_df = None


############################################################
################## DramV annotations  ######################
############################################################

print('Getting DRAMv annotations')
DRAMv_df = pd.DataFrame()
for sampleID in sample_list:
        sub_annotation = pd.read_csv(os.path.join(annotation_dir, 'DRAMv_annotate', sampleID + '.dramv_annotate.out', 'annotations.tsv'), sep='\t')
        sub_annotation['SequenceID'] = sampleID
        DRAMv_df = DRAMv_df.append(sub_annotation, ignore_index=True)
        sub_annotation = None

DRAMv_df.to_csv(os.path.join(annotation_dir, 'DRAMv_annotations.txt'), header=True, index=False, sep='\t')

DRAMv_df = None


############################################################
################### DramV distillates  #####################
############################################################

print('Getting DRAMv distillates')
DRAMv_df = pd.DataFrame()
for sampleID in sample_list:
        sub_annotation = pd.read_csv(os.path.join(annotation_dir, 'DRAMv_distill', sampleID + '.dramv_distilled.out', 'vMAG_stats.tsv'), sep='\t')
        sub_annotation['SequenceID'] = sampleID
        DRAMv_df = DRAMv_df.append(sub_annotation, ignore_index=True)
        sub_annotation = None

DRAMv_df.to_csv(os.path.join(annotation_dir, 'DRAMv_distallates.txt'), header=True, index=False, sep='\t')

DRAMv_df = None
