#!/usr/bin/env python

# score_var.py scores strcutural variant for disruption to 3D genome folding using Akita.
# Written in Python v 3.7.11

import akita_utils as utils
import numpy as np
import pandas as pd
import pysam
import argparse
import sys
import os
from pathlib import Path



# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Parse through arguments


parser = argparse.ArgumentParser(
                    prog = 'Akita_variant_scoring',
                    description='Score variants based on disruption to genome folding in the surround 1 Mb region using Akita (Fudenber et. al. 2020).')

parser.add_argument('--in',
                    dest = 'in_file',
                    help = 'Input file with variants.', 
                    type = str,
                    required = True)
parser.add_argument('--format',
                    dest = 'file_format', 
                    help = 'Format for input file.', 
                    type = str,
                    choices = ['vcf', 'df'],
                    default = 'vcf',
                    required = False)
parser.add_argument('--type',
                    dest = 'var_type', 
                    help = 'Variant type: simple or SV.', 
                    type = str,
                    choices = ['simple', 'SV'],
                    default = 'SV',
                    required = False)
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
                    help = 'Method(s) used to calculate disruption scores.', 
                    type = str,
                    choices = ['mse', 'corr', 'both'],
                    default = 'both',
                    required = False)
parser.add_argument('--shift_by', 
                    dest = 'shift_window',
                    nargs = '+', 
                    help = 'Values for shifting prediciton windows (e.g. to predict with no shift (variant centered) and shift by 1 bp to either side, use: -1 0 1). Values outside of range -450000 ≤ x ≤ 450000 will be ignored.',
                    type = int,
                    default = '0',
                    required = False)
parser.add_argument('--out',
                    dest = 'out_file', 
                    help = 'Directory in which to save the output file(s).', 
                    type = str,
                    default = 'output/out_file',
                    required = False)
parser.add_argument('--limit',
                    dest = 'svlen_limit', 
                    help = 'Maximum structural variant length to be scored.', 
                    type = int,
                    default = 700000,
                    required = False)

args = parser.parse_args()


in_file = args.in_file
file_format = args.file_format
var_type = args.var_type
fasta = args.fasta
chrom_lengths = args.chrom_lengths
centromere_coords = args.centromere_coords
score = args.score
shift_by = args.shift_window
out_file = args.out_file
svlen_limit = args.svlen_limit




# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Get necessary files if they are not there


if chrom_lengths == 'data/chrom_lengths_hg38' and not Path(chrom_lengths).is_file():
    os.system('wget -P ./data/ https://raw.githubusercontent.com/ketringjoni/Akita_variant_scoring/main/data/chrom_lengths_hg38')
        
if centromere_coords == 'data/centromere_coords_hg38' and not Path(centromere_coords).is_file():
    os.system('wget -P ./data/ https://raw.githubusercontent.com/ketringjoni/Akita_variant_scoring/main/data/centromere_coords_hg38')

if centromere_coords == 'data/hg38.fa' and not Path(fasta).is_file():
    os.system('wget -P ./data/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')
    os.system('gunzip data/hg38.fa.gz')

    
    
outdir = 'score_var_output'
if not os.path.exists(outdir):
    os.mkdir(outdir)
out_file = os.path.join(outdir, out_file)
    
    

    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Read in (and adjust) data

    
# Read in variants
variants = utils.read_input(in_file, file_format, var_type)

# Index input based on row number and create output with same indexes
variants['var_index'] = variants.index
variant_scores = pd.DataFrame({'var_index':variants.var_index})


# Read other necessary data
fasta_open = pysam.Fastafile(fasta) # Note: if you want to get sequence starting at POS you should use this command on POS-1 since the genome starts at 0 with this command but 1 with how positions are annotated
chrom_lengths = pd.read_table(chrom_lengths, header = None, names = ['CHROM', 'chrom_max'])
centromere_coords = pd.read_table(centromere_coords, header = None, names = ['CHROM', 'centro_start', 'centro_stop'])
nt = ['A', 'T', 'C', 'G']
max_shift = 450000


# Adjust shift input
if shift_by == 0:
    shift_by = [0]
else:
    # Remove shifts that are outside of allowed range
    shift_by = [x for x in shift_by if x > -max_shift and x < max_shift]

    
# Adjust SV length limit
if svlen_limit > 700000:
    svlen_limit = 700000
    

    


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Filter out variants that cannot be scored

# Get indexes for variants to exclude

# Exclude mitochondrial variants
chrM_var = pd.DataFrame({'var_index' : list(variants[variants.CHROM == 'chrM'].var_index),
                         'reason' : ' Mitochondrial chromosome.'})

# Exclude variants larger than limit
if var_type == 'SV':
    too_log_var = pd.DataFrame({'var_index' : [y for x,y in zip(variants.SVLEN, variants.var_index) 
                                 if not pd.isnull(x) and abs(int(x)) > svlen_limit],
                                'reason' : f' SV longer than {svlen_limit}.'})
else:
    too_log_var = pd.DataFrame()

filtered_out = pd.concat([chrM_var, too_log_var], axis = 0)
filtered_out.var_index = filtered_out.var_index.astype('int')

# Save filtered out variants into file
filtered_out.to_csv(f'{out_file}_filtered_out', sep = ':', index = False, header = False)

# Exclude
variants = variants[[x not in filtered_out.var_index.values for x in variants.var_index]]





# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Run: Make Akita predictions and calculate disruption scores


# Create log file to save standard output with error messages
std_output = sys.stdout
log_file = open(f'{out_file}_std_output.log','w')
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

            if var_type == 'simple' or all([x in nt for x in ALT]):
                scores = utils.get_scores(CHR, POS, REF, ALT, score, shift, chrom_lengths, centromere_coords, fasta_open)

            elif var_type == 'SV':
                END = variant.END
                SVTYPE = variant.SVTYPE

                scores = utils.get_scores_SV(CHR, POS, ALT, END, SVTYPE, score, shift, 
                                             chrom_lengths, centromere_coords, fasta_open) 


            if score == 'corr':
                variant_scores.loc[variant_scores.var_index == var_index, f'corr_{shift}'] = scores
            elif score == 'mse':
                variant_scores.loc[variant_scores.var_index == var_index, f'MSE_{shift}'] = scores
            elif score == 'both':
                variant_scores.loc[variant_scores.var_index == var_index, f'MSE_{shift}'] = scores[0] 
                variant_scores.loc[variant_scores.var_index == var_index, f'corr_{shift}'] = scores[1]

            print(str(var_index) + ' (' + str(shift) + ' shift)')

        except Exception as e: 

            print(str(var_index) + ' (' + str(shift) + ' shift)' + ':', e)

            pass



# Save results
variant_scores.to_csv(f'{out_file}_scores', sep = '\t', index = False)


# Save standard output with error messages to log file
sys.stdout = std_output
log_file.close()






