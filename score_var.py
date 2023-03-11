#!/usr/bin/env python

# score_var.py scores strcutural variant for disruption to 3D genome folding using Akita.
# Written in Python v 3.7.11

# Use score_var.py <variant_df> <fasta> <chrom_lengths> <centromere_coords> <scores> <out_dir>

import akita_utils as utils
import numpy as np
import pandas as pd
import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("variant_df", help="Dataframe generated from vcf file with variant info.", type=str)
parser.add_argument("fasta", help="hg38 reference genome fasta file.", type=str)
parser.add_argument("chrom_lengths", help="File with lengths of chromosomes in hg38. Columns: chromosome (ex: 1), length; no header.", type=str)
parser.add_argument("centromere_coords", help="Centromere coordinates for hg38. Can get here: https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz. Columns: chromosome (ex: chr1), start, end; no header.", type=str)
parser.add_argument("scores", help="Method used to calculate disruption scores. Default: both", 
    choices=['MSE', 'correlation', 'both'])
parser.add_argument("out_dir", help="Directory in which to save the output file(s).", type=str)

args = parser.parse_args()


variant_df = args.variant_df
fasta = args.fasta
chrom_lengths = args.chrom_lengths
centromere_coords = args.centromere_coords
scores = args.scores
out_dir = args.out_dir


# Read variants, fasta, chrom_lengths, centromere_coords

variants = pd.read_csv(variant_df, sep = '\t')
variants = variants[variants.CHROM != 'chrM'] # Exclude mitochondrial variants

fasta_open = pysam.Fastafile(fasta) # Note: if you want to get sequence starting at POS you should use this command on POS-1 since the genome starts at 0 with this command but 1 with how positions are annotated

chrom_lengths = pd.read_table(chrom_lengths, header = None, names = ['CHROM', 'chrom_max'])

centromere_coords = pd.read_table(centromere_coords, header = None, names = ['CHROM', 'centro_start', 'centro_stop'])

nt = ['A', 'T', 'C', 'G']


# Make Akita predictions and calculate disruption scores
for i in variants.index: 
    
    variant = variants.loc[i]

    CHR = variant.CHROM
    POS = variant.POS
    REF = variant.REF
    ALT = variant.ALT
    if 'SVTYPE' in variant.index:
        SVTYPE = variant.SVTYPE
    
    print(str(i+1)+'/'+str(len(variants)))
    
    # Get disruption scores

    # Simple variants
    if all([x in nt for x in ALT]):
        
        try:
            if scores == 'both':
                MSE, correlation = utils.get_scores(CHR, POS, REF, ALT, chrom_lengths, centromere_coords, fasta_open)
        except:
            print(str(i)+': Variant cannot be scored likely due to length or N composition')
            pass
          
    # Structural variants
    elif SVTYPE in ['DEL', 'DUP', 'INV', 'BND']:
        if SVTYPE != 'BND':
            END = int(variant.END)
        else:
            END = np.nan

        try:
            if scores == 'both':
                MSE, correlation = utils.get_scores_SV(CHR, POS, ALT, END, SVTYPE, chrom_lengths, centromere_coords, fasta_open) 
        except:
            print(str(i)+': Variant cannot be scored likely due to length or N composition')
            pass
        
    else:
        print(str(i)+': SV type or input format not supported. Note: Predictions cannot be made for insertions nor non-tandem duplications.')
        MSE, correlation = np.nan, np.nan
         
            
    variants.loc[i, 'MSE'] = MSE
    variants.loc[i, 'corr'] = correlation


# Save results
variants.to_csv(out_dir, sep = '\t', index = False)





