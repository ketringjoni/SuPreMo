#!/usr/bin/env python

# score_var.py scores strcutural variant for disruption to 3D genome folding using Akita.
# Written in Python v 3.7.11

import akita_utils as utils
import numpy as np
import pandas as pd
import argparse

nt = ['A', 'T', 'C', 'G']

parser = argparse.ArgumentParser()

parser.add_argument("variant_df", help="Dataframe generated from vcf file with variant info.", type=str)
parser.add_argument("scores", help="Method used to calculate disruption scores. Default: both", 
    choices=['MSE', 'correlation', 'both'])
parser.add_argument("out_dir", help="Directory in which to save the output file(s).", type=str)
args = parser.parse_args()


variant_df = args.variant_df
scores = args.scores
out_dir = args.out_dir


# Read variants 
variants = pd.read_csv(variant_df, sep = '\t')

# Exclude mitochondrial variants
variants = variants[variants.CHROM != 'chrM']


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
    try:
          # Simple variants
        if all([x in nt for x in ALT]):
              if scores == 'both':
                    MSE, correlation = utils.get_scores(CHR, POS, REF, ALT)
          
          
          # Structural variants
        elif SVTYPE in ['INS']:
            raise ValueError('Predictions cannot be made for insertions nor non-tandem duplications.')
          
        elif SVTYPE in ['DEL', 'DUP', 'INV', 'BND']:
            if SVTYPE != 'BND':
                END = int(variant.END)
            else:
                END = np.nan
                    
            if scores == 'both':
                  MSE, correlation = utils.get_scores_SV(CHR, POS, ALT, END, SVTYPE) 
        else:
            raise ValueError('SV type or input format not supported')
            
        variants.loc[i, 'MSE'] = MSE
        variants.loc[i, 'corr'] = correlation

        
    except:
        pass


# Save results
variants.to_csv(out_dir, sep = '\t', index = False)





