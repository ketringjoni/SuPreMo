#!/usr/bin/env python

# Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020).
# Written in Python v 3.7.11

'''
usage: Akita_variant_scoring [-h] --in IN_FILE [--fa FASTA]
[--chrom CHROM_LENGTHS]
                             [--centro CENTROMERE_COORDS]
                             [--score {mse,corr} [{mse,corr} ...]]
                             [--shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]]
                             [--file OUT_FILE] [--dir OUT_DIR]
                             [--limit SVLEN_LIMIT]

Pipeline for scoring variants for disruption to genome folding using Akita
(Fudenberg et. al. 2020).

optional arguments:
  -h, --help            show this help message and exit
  --in IN_FILE          Input file with variants.
  --fa FASTA            hg38 reference genome fasta file.
  --chrom CHROM_LENGTHS
                        File with lengths of chromosomes in hg38. Columns:
                        chromosome (ex: 1), length; no header.
  --centro CENTROMERE_COORDS
                        Centromere coordinates for hg38. Columns: chromosome
                        (ex: chr1), start, end; no header.
  --score {mse,corr} [{mse,corr} ...]
                        Method(s) used to calculate disruption scores.
  --shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]
                        Values for shifting prediciton windows (e.g. to
                        predict with no shift (variant centered) and shift by
                        1 bp to either side, use: -1 0 1). Values outside of
                        range -450000 ≤ x ≤ 450000 will be ignored. Prediction
                        windows at the edge of chromosome arms will only be
                        shifted in the direction that is possible (ex. for
                        window at chrom start, a -1 shift will be treated as a
                        1 shift since it is not possible to shift left.)
  --file OUT_FILE       Prefix for output files.
  --dir OUT_DIR         Output directory.
  --limit SVLEN_LIMIT   Maximum structural variant length to be scored (<=
                        700000).
                        
'''



# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Parse through arguments

import argparse

parser = argparse.ArgumentParser(
                    prog = 'Akita_variant_scoring',
                    description='Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020).')

parser.add_argument('--in',
                    dest = 'in_file',
                    help = 'Input file with variants.', 
                    type = str,
                    required = True)

parser.add_argument('--fa',
                    dest = 'fasta', 
                    help = 'hg38 reference genome fasta file.', 
                    type = str,
                    default = 'data/hg38.fa',
                    required = False)

parser.add_argument('--chrom',
                    dest = 'chrom_lengths', 
                    help = 'File with lengths of chromosomes in hg38. Columns: chromosome (ex: 1), length; no header.', 
                    type = str,
                    default = 'data/chrom_lengths_hg38',
                    required = False)

parser.add_argument('--centro',
                    dest = 'centromere_coords', 
                    help = 'Centromere coordinates for hg38. Columns: chromosome (ex: chr1), start, end; no header.', 
                    type = str,
                    default = 'data/centromere_coords_hg38',
                    required = False)

parser.add_argument('--score',
                    dest = 'score', 
                    nargs = '+', 
                    help = 'Method(s) used to calculate disruption scores.', 
                    type = str,
                    choices = ['mse', 'corr'],
                    default = ['mse', 'corr'],
                    required = False)

parser.add_argument('--shift_by', 
                    dest = 'shift_window',
                    nargs = '+', 
                    help = 'Values for shifting prediciton windows (e.g. to predict with no shift (variant centered) and shift by 1 bp to either side, use: -1 0 1). Values outside of range -450000 ≤ x ≤ 450000 will be ignored. Prediction windows at the edge of chromosome arms will only be shifted in the direction that is possible (ex. for window at chrom start, a -1 shift will be treated as a 1 shift since it is not possible to shift left.)',
                    type = int,
                    default = [0],
                    required = False)

parser.add_argument('--file',
                    dest = 'out_file', 
                    help = 'Prefix for output files.', 
                    type = str,
                    default = 'score_var_results',
                    required = False)
parser.add_argument('--dir',
                    dest = 'out_dir', 
                    help = 'Output directory.', 
                    type = str,
                    default = 'score_var_output',
                    required = False)

def max_svlen_limit(x):
    x = int(x)
    if x > 700000:
        raise argparse.ArgumentTypeError("Maximum SV length limit is 700000.")
    return x

parser.add_argument('--limit',
                    dest = 'svlen_limit', 
                    help = 'Maximum structural variant length to be scored (<= 700000).', 
                    type = max_svlen_limit,
                    default = 700000,
                    required = False)

args = parser.parse_args()


in_file = args.in_file
fasta_path = args.fasta
chrom_lengths_path = args.chrom_lengths
centromere_coords_path = args.centromere_coords
scores = args.score
shift_by = args.shift_window
out_file = args.out_file
out_dir = args.out_dir
svlen_limit = args.svlen_limit





# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Get necessary files if they are not there

import os
from pathlib import Path


if chrom_lengths_path == 'data/chrom_lengths_hg38' and not Path(chrom_lengths_path).is_file():
    os.system('wget -P ./data/ https://raw.githubusercontent.com/ketringjoni/Akita_variant_scoring/main/data/chrom_lengths_hg38')
    print('Chromosome lengths file downloaded as data/chrom_lengths_hg38.')
        
if centromere_coords_path == 'data/centromere_coords_hg38' and not Path(centromere_coords_path).is_file():
    os.system('wget -P ./data/ https://raw.githubusercontent.com/ketringjoni/Akita_variant_scoring/main/data/centromere_coords_hg38')
    print('Centromere coordinates file downloaded as data/centromere_coords_hg38.')

if fasta_path == 'data/hg38.fa' and not Path(fasta_path).is_file():
    os.system('wget -P ./data/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')
    os.system('gunzip data/hg38.fa.gz')
    print('Fasta file downloaded as data/hg38.fa.')

    
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
out_file = os.path.join(out_dir, out_file)
 
    
    

    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Read in (and adjust) data

import akita_utils as utils
import pandas as pd
import pysam
   
# Read in variants
variants = utils.read_input(in_file)

# Index input based on row number and create output with same indexes
variants['var_index'] = variants.index
variant_scores = pd.DataFrame({'var_index':variants.var_index})


# Read other necessary data
fasta_open = pysam.Fastafile(fasta_path)
chrom_lengths = pd.read_table(chrom_lengths_path, header = None, names = ['CHROM', 'chrom_max'])
centromere_coords = pd.read_table(centromere_coords_path, header = None, names = ['CHROM', 'centro_start', 'centro_stop'])

utils.fasta_open = fasta_open
utils.chrom_lengths = chrom_lengths
utils.centromere_coords = centromere_coords



# Adjust shift input: Remove shifts that are outside of allowed range
shift_by = [x for x in shift_by if x > -utils.max_shift and x < utils.max_shift]
    
  

    


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Filter out variants that cannot be scored

# Get indexes for variants to exclude

# Exclude mitochondrial variants
chrM_var = pd.DataFrame({'var_index' : list(variants[variants.CHROM == 'chrM'].var_index),
                         'reason' : ' Mitochondrial chromosome.'})

# Exclude variants larger than limit
if 'SVLEN' in variants.columns:
    too_long_var = pd.DataFrame({'var_index' : [y for x,y in zip(variants.SVLEN, variants.var_index) 
                                                if not pd.isnull(x) and abs(int(x)) > svlen_limit],
                                 'reason' : f' SV longer than {svlen_limit}.'})
else:
    too_long_var = pd.DataFrame()

filtered_out = pd.concat([chrM_var, too_long_var], axis = 0)
filtered_out.var_index = filtered_out.var_index.astype('int')

# Save filtered out variants into file
filtered_out.to_csv(f'{out_file}_filtered_out', sep = ':', index = False, header = False)

# Exclude
variants = variants[[x not in filtered_out.var_index.values for x in variants.var_index]]





# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Run: Make Akita predictions and calculate disruption scores

import sys
import numpy as np

# Create log file to save standard output with error messages
std_output = sys.stdout
log_file = open(f'{out_file}_log','w')
sys.stdout = log_file


# Loop through each row (not index) and get disruption scores 
for i in range(len(variants)):
    
    variant = variants.iloc[i]

    var_index = variant.var_index
    CHR = variant.CHROM
    POS = variant.POS
    REF = variant.REF
    ALT = variant.ALT

    for shift in shift_by:
        
        try:

            if all([x in utils.nt for x in REF]) and all([x in utils.nt for x in ALT]):
                scores = utils.get_scores(CHR, POS, REF, ALT, scores, shift)

            else:
                END = variant.END
                SVTYPE = variant.SVTYPE

                scores = utils.get_scores_SV(CHR, POS, ALT, END, SVTYPE, scores, shift) 

            for score in scores:
                variant_scores.loc[variant_scores.var_index == var_index, f'{score}_{shift}'] = scores[score]

            print(str(var_index) + ' (' + str(shift) + ' shift)')

        except Exception as e: 

            print(str(var_index) + ' (' + str(shift) + ' shift)' + ': Error:', e)

            pass



# Save results
variant_scores.to_csv(f'{out_file}_scores', sep = '\t', index = False)


# Save standard output with error messages to log file
sys.stdout = std_output
log_file.close()


# Adjust log file to only have 1 row per variant
if os.path.exists(f'{out_file}_log'):

    log_file = pd.read_csv(f'{out_file}_log', sep = '\n', names = ['output']) 
    # Move warnings (printed 1 line before variant) to variant line
    indexes = np.array([[index, index+1] for (index, item) in enumerate(log_file.output) if item.startswith('Warning')])
    
    if len(indexes) != 0:
        log_file.loc[indexes[:,1],'output'] = [x+': '+y for x,y in zip(list(log_file.loc[indexes[:,1],'output']),
                                                                       list(log_file.loc[indexes[:,0],'output']))]
        log_file.drop(indexes[:,0], axis = 0, inplace = True)
        log_file.to_csv(f'{out_file}_log', sep = '\t', header = None, index = False)



