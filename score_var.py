#!/usr/bin/env python

# Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020).
# Written in Python v 3.7.11

'''
usage: Akita_variant_scoring [-h] --in IN_FILE [--fa FASTA]
                             [--scores {mse,corr,ssi,scc,ins,di,dec,tri,pca} [{mse,corr,ssi,scc,ins,di,dec,tri,pca} ...]]
                             [--shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]]
                             [--file OUT_FILE] [--dir OUT_DIR]
                             [--limit SVLEN_LIMIT] [--seq_len SEQ_LEN]
                             [--revcomp {no_revcomp,add_revcomp,only_revcomp}]
                             [--augment] [--get_seq] [--get_tracks]
                             [--get_maps] [--get_scores] [--nrows NROWS]

Pipeline for generating mutated sequences for input into predictive models and for scoring variants for disruption to genome folding.

optional arguments:
  -h, --help            show this help message and exit
  --in IN_FILE          
                        Input file with variants. Accepted formats are:  
                            VCF file,
                            TSV file (SVannot output), 
                            BED file,
                            tab-delimited text file.
                        Can be gzipped. Coordinates should be 1-based and left open (,], except for SNPs (follow vcf 4.1/4.2 specifications). See custom_perturbations.ipynb for BED or TXT file input specifications.
                         (default: None)
  --fa FASTA            Optional path to hg38 reference genome fasta file. If not provided and not existing in data/, it will be downloaded.
                         (default: data/hg38.fa)
  --scores {mse,corr,ssi,scc,ins,di,dec,tri,pca} [{mse,corr,ssi,scc,ins,di,dec,tri,pca} ...]
                        
                        Method(s) used to calculate disruption scores. Use abbreviations as follows:
                            mse: Mean squared error
                            corr: Spearman correlation
                            ssi: Structural similarity index measure
                            scc: Stratum adjusted correlation coefficient
                            ins: Insulation
                            di: Directionality index
                            dec: Contact decay
                            tri: Triangle method
                            pca: Principal component method.
                         (default: ['mse', 'corr'])
  --shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]
                        Values for shifting prediciton windows inputted as space-separated integers (e.g. -1 0 1). Values outside of range -450000 ≤ x ≤ 450000 will be ignored. Prediction windows at the edge of chromosome arms will only be shifted in the direction that is possible (ex. for window at chrom start, a -1 shift will be treated as a 1 shift since it is not possible to shift left).
                         (default: [0])
  --file OUT_FILE       Prefix for output files. Saved files will overwrite any existing files.
                         (default: score_var_results)
  --dir OUT_DIR         Output directory. (default: score_var_output)
  --limit SVLEN_LIMIT   Maximum length of variants to be scored. Filtering out variants that are too big can save time and memory. If not specified, will be set to 2/3 of seq_len.
                         (default: None)
  --seq_len SEQ_LEN     Length for sequences to generate. Default value is based on Akita requirement. If non-default value is set, get_scores must be false.
                         (default: 1048576)
  --revcomp {no_revcomp,add_revcomp,only_revcomp}
                        
                        Option to use the reverse complement of the sequence:
                            no_revcomp: no, only use the standard sequence;
                            add_revcomp: yes, use both the standard sequence and its reverse complement;
                            only_revcomp: yes, only use the reverse complement of the sequence.
                        
                         (default: ['no_revcomp'])
  --augment             
                        Only applicable if --get_scores is specified. Get the average score from four sequences: 
                            1) no augmentation: 0 shift and no reverse complement, 
                            2) +1bp shift and no reverse complement, 
                            3) -1bp shift and no reverse complement, 
                            4) 0 shift and take reverse complement. 
                            Overwrites --shift_by and --revcomp arguments.
                         (default: False)
  --get_seq             Save sequences for the reference and alternate alleles in fa file format. If --get_seq is not specified, must specify --get_scores.
                        
                        Sequence name format: {var_index}_{shift}_{revcomp_annot}_{seq_index}_{var_rel_pos}
                            var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles; 
                            shift: integer that window is shifted by; 
                            revcomp_annot: present only if reverse complement of sequence was taken; 
                            seq_index: index for sequences generated for that variant: 0-1 for non-BND reference and alternate sequences and 0-2 for BND left and right reference sequence and alternate sequence; 
                            var_rel_pos: relative position of variant in sequence: list of two for non-BND variant positions in reference and alternate sequence and an integer for BND breakend position in reference and alternate sequences. 
                        
                        There are 2-3 entries per prediction (2 for non-BND variants and 3 for BND variants).
                        
                        To read fasta file: pysam.Fastafile(filename).fetch(seqname, start, end).upper(). 
                        To get sequence names in fasta file: pysam.Fastafile(filename).references.
                         (default: False)
  --get_tracks          Save disruption score tracks (448 bins) in npy file format. Only possible for mse and corr scores. --get_scores must be specified.
                                            
                        Dictionary item name format: {var_index}_{track}_{shift}_{revcomp_annot}
                            var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles; 
                            track: disruption score track specified; 
                            shift: integer that window is shifted by; 
                            revcomp_annot: present only if reverse complement of sequence was taken. 
                            
                        There is 1 entry per prediction: a 448x1 array.
                                            
                        To read into a dictionary in python: np.load(filename, allow_pickle="TRUE").item()
                         (default: False)
  --get_maps            Save predicted contact frequency maps in npy file format. --get_scores must be specified.
                                            
                        Dictionary item name format: {var_index}_{shift}_{revcomp_annot}
                            var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles; 
                            shift: integer that window is shifted by; 
                            revcomp_annot: present only if reverse complement of sequence was taken. 
                                            
                        There is 1 entry per prediction. Each entry has 2-3 448x448 arrays (2 for non-BND variants and 3 for BND variants). 
                        
                        To read into a dictionary in python: np.load(filename, allow_pickle="TRUE").item()
                         (default: False)
  --get_scores          Get disruption scores. If --get_scores is not specified, must specify --get_seq. Scores saved in a dataframe with the same number of rows as the input. For multiple alternate alleles, the scores are separated by a comma. To convert the scores from strings to integers, use float(x), after separating rows with multiple alternate alleles. Scores go up to 20 decimal points.
                         (default: False)
  --nrows NROWS         Number of rows (perturbations) to read at a time from input. When dealing with large inputs, selecting a subset of rows to read at a time allows scores to be saved in increments and uses less memory. Files with scores and filtered out variants will be temporarily saved in output direcotry. The file names will have a suffix corresponding to the set of nrows (0-based), for example for an input with 2700 rows and with nrows = 1000, there will be 3 sets. At the end of the run, these files will be concatenated into a comprehensive file and the temporary files will be removed.
                                             (default: 1000)
                        
'''



# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Parse through arguments

import argparse

class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

parser = argparse.ArgumentParser(
                    prog = 'Akita_variant_scoring',
                    formatter_class = CustomFormatter,
                    description='''Pipeline for generating mutated sequences for input into predictive models and for scoring variants for disruption to genome folding.''')

parser.add_argument('--in',
                    dest = 'in_file',
                    help = '''
Input file with variants. Accepted formats are:  
    VCF file,
    TSV file (SVannot output), 
    BED file,
    tab-delimited text file.
Can be gzipped. Coordinates should be 1-based and left open (,], except for SNPs (follow vcf 4.1/4.2 specifications). See custom_perturbations.ipynb for BED or TXT file input specifications.
''', 
                    type = str,
                    required = True)

parser.add_argument('--fa',
                    dest = 'fasta', 
                    help = '''Optional path to hg38 reference genome fasta file. If not provided and not existing in data/, it will be downloaded.
''', 
                    type = str,
                    default = 'data/hg38.fa',
                    required = False)

parser.add_argument('--scores',
                    dest = 'scores', 
                    nargs = '+', 
                    help = '''
Method(s) used to calculate disruption scores. Use abbreviations as follows:
    mse: Mean squared error
    corr: Spearman correlation
    ssi: Structural similarity index measure
    scc: Stratum adjusted correlation coefficient
    ins: Insulation
    di: Directionality index
    dec: Contact decay
    tri: Triangle method
    pca: Principal component method.
''', 
                    type = str,
                    choices = ['mse', 'corr',  'ssi', 'scc', 'ins',  'di', 'dec', 'tri', 'pca'],
                    default = ['mse', 'corr'],
                    required = False)

parser.add_argument('--shift_by', 
                    dest = 'shift_window',
                    nargs = '+', 
                    help = '''Values for shifting prediciton windows inputted as space-separated integers (e.g. -1 0 1). Values outside of range -450000 ≤ x ≤ 450000 will be ignored. Prediction windows at the edge of chromosome arms will only be shifted in the direction that is possible (ex. for window at chrom start, a -1 shift will be treated as a 1 shift since it is not possible to shift left).
''',
                    type = int,
                    default = [0],
                    required = False)

parser.add_argument('--file',
                    dest = 'out_file', 
                    help = '''Prefix for output files. Saved files will overwrite any existing files.
''', 
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
                    help = '''Maximum length of variants to be scored. Filtering out variants that are too big can save time and memory. If not specified, will be set to 2/3 of seq_len.
''', 
                    type = int,
                    required = False)

parser.add_argument('--seq_len',
                    dest = 'seq_len', 
                    help = '''Length for sequences to generate. Default value is based on Akita requirement. If non-default value is set, get_scores must be false.
''', 
                    type = int,
                    default = 1048576,
                    required = False)

parser.add_argument('--revcomp',
                    dest = 'revcomp',
                    nargs = 1, 
                    help = '''
Option to use the reverse complement of the sequence:
    no_revcomp: no, only use the standard sequence;
    add_revcomp: yes, use both the standard sequence and its reverse complement;
    only_revcomp: yes, only use the reverse complement of the sequence.

''', 
                    type = str,
                    choices = ['no_revcomp', 'add_revcomp', 'only_revcomp'],
                    default = ['no_revcomp'],
                    required = False)

parser.add_argument('--augment',
                    dest = 'augment', 
                    help = '''
Only applicable if --get_scores is specified. Get the average score from four sequences: 
    1) no augmentation: 0 shift and no reverse complement, 
    2) +1bp shift and no reverse complement, 
    3) -1bp shift and no reverse complement, 
    4) 0 shift and take reverse complement. 
    Overwrites --shift_by and --revcomp arguments.
''', 
                    action='store_true',
                    required = False)

parser.add_argument('--get_seq',
                    dest = 'get_seq', 
                    help = '''Save sequences for the reference and alternate alleles in fa file format. If --get_seq is not specified, must specify --get_scores.

Sequence name format: {var_index}_{shift}_{revcomp_annot}_{seq_index}_{var_rel_pos}
    var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles; 
    shift: integer that window is shifted by; 
    revcomp_annot: present only if reverse complement of sequence was taken; 
    seq_index: index for sequences generated for that variant: 0-1 for non-BND reference and alternate sequences and 0-2 for BND left and right reference sequence and alternate sequence; 
    var_rel_pos: relative position of variant in sequence: list of two for non-BND variant positions in reference and alternate sequence and an integer for BND breakend position in reference and alternate sequences. 

There are 2-3 entries per prediction (2 for non-BND variants and 3 for BND variants).

To read fasta file: pysam.Fastafile(filename).fetch(seqname, start, end).upper(). 
To get sequence names in fasta file: pysam.Fastafile(filename).references.
''', 
                    action = 'store_true',
                    required = False)

parser.add_argument('--get_tracks',
                    dest = 'get_tracks', 
                    help = '''Save disruption score tracks (448 bins) in npy file format. Only possible for mse and corr scores. --get_scores must be specified.
                    
Dictionary item name format: {var_index}_{track}_{shift}_{revcomp_annot}
    var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles; 
    track: disruption score track specified; 
    shift: integer that window is shifted by; 
    revcomp_annot: present only if reverse complement of sequence was taken. 
    
There is 1 entry per prediction: a 448x1 array.
                    
To read into a dictionary in python: np.load(filename, allow_pickle="TRUE").item()
''', 
                    action='store_true', 
                    required = False)

parser.add_argument('--get_maps',
                    dest = 'get_maps', 
                    help = '''Save predicted contact frequency maps in npy file format. --get_scores must be specified.
                    
Dictionary item name format: {var_index}_{shift}_{revcomp_annot}
    var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles; 
    shift: integer that window is shifted by; 
    revcomp_annot: present only if reverse complement of sequence was taken. 
                    
There is 1 entry per prediction. Each entry has 2-3 448x448 arrays (2 for non-BND variants and 3 for BND variants). 

To read into a dictionary in python: np.load(filename, allow_pickle="TRUE").item()
''',
                    action='store_true', 
                    required = False)

parser.add_argument('--get_scores',
                    dest = 'get_scores', 
                    help = '''Get disruption scores. If --get_scores is not specified, must specify --get_seq. Scores saved in a dataframe with the same number of rows as the input. For multiple alternate alleles, the scores are separated by a comma. To convert the scores from strings to integers, use float(x), after separating rows with multiple alternate alleles. Scores go up to 20 decimal points.
''', 
                    action='store_true',
                    required = False)

parser.add_argument('--nrows',
                    dest = 'nrows', 
                    help = '''Number of rows (perturbations) to read at a time from input. When dealing with large inputs, selecting a subset of rows to read at a time allows scores to be saved in increments and uses less memory. Files with scores and filtered out variants will be temporarily saved in output direcotry. The file names will have a suffix corresponding to the set of nrows (0-based), for example for an input with 2700 rows and with nrows = 1000, there will be 3 sets. At the end of the run, these files will be concatenated into a comprehensive file and the temporary files will be removed.
                    ''', 
                    type = int,
                    default = 1000,
                    required = False)


args = parser.parse_args()


in_file = args.in_file
fasta_path = args.fasta
scores_to_use = args.scores
shift_by = args.shift_window
out_file = args.out_file
out_dir = args.out_dir
svlen_limit = args.svlen_limit
seq_len = args.seq_len
revcomp = args.revcomp[0]
augment = args.augment
get_seq = args.get_seq
get_tracks = args.get_tracks
get_maps = args.get_maps
get_scores = args.get_scores
var_set_size = args.nrows


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Adjust inputs from arguments


# Handle argument dependencies

if seq_len != 1048576:
    get_scores = False

if svlen_limit is None:
    svlen_limit = 2/3*seq_len

if not get_seq and not get_scores:
    raise ValueError('Either get_seq and/or get_scores must be True.')
    

# Adjust shift input: Remove shifts that are outside of allowed range
max_shift = 0.4*seq_len
shift_by = [x for x in shift_by if x > -max_shift and x < max_shift]


# Adjust input for taking the reverse complement
if revcomp == 'no_revcomp':
    revcomp_decision = [False]
elif revcomp == 'add_revcomp':
    revcomp_decision = [False, True]
elif revcomp == 'only_revcomp':
    revcomp_decision = [True]

    
# Adjust input for taking the average score from augmented sequences
if augment:
    shift_by = [-1,0,1]
    revcomp_decision = [False, True]
#     scores_to_use = [x for x in scores_to_use if x in ['mse', 'corr']]

revcomp_decision_i = revcomp_decision


# Create dictionaries to save sequences, maps, and disruption score tracks, if specified
if get_seq:
    sequences = {}
    
if get_maps:
    variant_maps = {}
    if get_scores == False:
        get_scores = True
        print('Must get scores to get maps. --get_scores was not specified but will be applied.')
    
if get_tracks:
    variant_tracks = {}
    if get_scores == False:
        get_scores = True
        print('Must get scores to get tracks. --get_scores was not specified but will be applied.')


    




# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Get necessary files if they are not there

import os
from pathlib import Path

chrom_lengths_path = 'data/chrom_lengths_hg38'
if not Path(chrom_lengths_path).is_file():
    os.system('wget -P ./data/ https://raw.githubusercontent.com/ketringjoni/Akita_variant_scoring/main/data/chrom_lengths_hg38')
    print('Chromosome lengths file downloaded as data/chrom_lengths_hg38.')

centromere_coords_path = 'data/centromere_coords_hg38'
if not Path(centromere_coords_path).is_file():
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

import pandas as pd
import pysam

chrom_lengths = pd.read_table(chrom_lengths_path, header = None, names = ['CHROM', 'chrom_max'])
centromere_coords = pd.read_table(centromere_coords_path, sep = '\t')
fasta_open = pysam.Fastafile(fasta_path)


# Assign necessary values to variables across module

# Module 1: reading utilities
import reading_utils
reading_utils.var_set_size = var_set_size


# Module 2: get_seq utilities
import get_seq_utils
get_seq_utils.fasta_open = fasta_open
get_seq_utils.chrom_lengths = chrom_lengths
get_seq_utils.centromere_coords = centromere_coords

get_seq_utils.svlen_limit = svlen_limit
get_seq_utils.seq_length = seq_len
get_seq_utils.half_patch_size = round(seq_len/2)


# Module 2: get_scores utilities
if get_scores:
    import get_scores_utils
    get_scores_utils.chrom_lengths = chrom_lengths
    get_scores_utils.centromere_coords = centromere_coords

    

    
    
    
import sys
import numpy as np

# Create log file to save standard output with error messages
print(f'Log file being saved here: {out_file}_log')
std_output = sys.stdout
log_file = open(f'{out_file}_log','w')
sys.stdout = log_file



    
var_set = 0
var_set_list = []
while True:
    
    # Read in variants
    variants = reading_utils.read_input(in_file, var_set)
    if len(variants) == 0:
        break
        
    # Index input based on row number and create output with same indexes
    variants['var_index'] = list(range(var_set*var_set_size, var_set*var_set_size + len(variants)))
    variants['var_index'] = variants['var_index'].astype(str)
    
    # If there are multiple alternate alleles, split those into new rows and indexes
    if any([',' in x for x in variants.ALT]):

        variants = (variants
                 .set_index(['CHROM', 'POS', 'REF', 'var_index'])
                 .apply(lambda x: x.str.split(',').explode())
                 .reset_index())

        g = variants.groupby(['var_index'])
        variants.loc[g['var_index'].transform('size').gt(1),
               'var_index'] += '-'+g.cumcount().astype(str)
        
    variant_scores = pd.DataFrame({'var_index':variants.var_index})
    
  


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
    filtered_out.var_index = filtered_out.var_index.astype('str')

    # Save filtered out variants into file
    filtered_out.to_csv(f'{out_file}_filtered_out_{var_set}', sep = ':', index = False, header = False)

    # Exclude
    variants = variants[[x not in filtered_out.var_index.values for x in variants.var_index]]



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Run: Make Akita predictions and calculate disruption scores



    # Loop through each row (not index) and get disruption scores 
    for i in range(len(variants)):

        variant = variants.iloc[i]

        var_index = variant.var_index
        CHR = variant.CHROM
        POS = variant.POS
        REF = variant.REF
        ALT = variant.ALT
            

        if 'SVTYPE' in variants.columns:
            END = int(variant.END)
            SVTYPE = variant.SVTYPE
        else:
            END = np.nan
            SVTYPE = np.nan

        for shift in shift_by:

            if augment: # if getting augmented score, take reverse complement only with 0 shift
                if shift != 0 & True in revcomp_decision:
                    revcomp_decision_i = [False]
                else:
                    revcomp_decision_i = revcomp_decision

            for revcomp in revcomp_decision_i:

                try:

                    if revcomp:
                        revcomp_annot = '_revcomp'
                    else:
                        revcomp_annot = ''
                    
                    sequences_i = get_seq_utils.get_sequences_SV(CHR, POS, REF, ALT, END, SVTYPE, shift, revcomp)


                    if get_seq:

                        # Get relative position of variant in sequence
                        var_rel_pos = str(sequences_i[-1]).replace(', ', '_')

                        for ii in range(len(sequences_i[:-1][:3])): 
                            sequences[f'{var_index}_{shift}{revcomp_annot}_{ii}_{var_rel_pos}'] = sequences_i[:-1][ii]

                    if get_scores:

                        scores = get_scores_utils.get_scores(CHR, POS, REF, ALT, sequences_i, SVTYPE, 
                                                  scores_to_use, shift, revcomp, get_tracks, get_maps)

                        if get_tracks:
                            for track in [x for x in scores.keys() if 'track' in x]:
                                variant_tracks[f'{var_index}_{track}_{shift}{revcomp_annot}'] = scores[track]
                                del scores[track]

                        if get_maps:
                            variant_maps[f'{var_index}_{shift}{revcomp_annot}'] = scores['maps']
                            del scores['maps']

                        for score in scores:
                            variant_scores.loc[variant_scores.var_index == var_index, 
                                               f'{score}_{shift}{revcomp_annot}'] = scores[score]


                    print(str(var_index) + ' (' + str(shift) + f' shift{revcomp_annot})')

                except Exception as e: 

                    print(str(var_index) + ' (' + str(shift) + f' shift{revcomp_annot})' + ': Error:', e)

                    pass
 
    
      
    # Combine results from all sets


    # Write sequences to fasta file
    if get_seq:

        if var_set == 0:
            sequences_all = sequences.copy()
        else:
            sequences_all.update(sequences)

    # Write scores to data frame
    if get_scores:

        # Take average of augmented sequences
        if augment:
            for score in scores:
                cols = [x for x in variant_scores.columns if score in x ]
                variant_scores[f'{score}_mean'] = variant_scores[cols].mean(axis = 1)
                variant_scores.drop(cols, axis = 1, inplace = True)

        # Convert scores from float to string so you can merge scores for variants with multiple alleles
        for col in variant_scores.iloc[:,1:].columns:
            variant_scores[col] = [format(x, '.20f') for x in variant_scores[col]]
            
                
        # Join scores for alternate alleles, separated by a comma
        if any(['-' in x for x in variant_scores.var_index]):
            variant_scores['var_index'] = variant_scores.var_index.str.split('-').str[0]
            
            variant_scores = (variant_scores
                              .set_index(['var_index'], drop = False)
                              .rename(columns = {'var_index':'var_index2'})
                              .groupby('var_index2')
                              .transform(','.join)
                              .reset_index()
                              .drop_duplicates())
        
        if var_set == 0:
            variant_scores.to_csv(f'{out_file}_scores_{var_set}', sep = '\t', index = False)
        else:
            variant_scores.to_csv(f'{out_file}_scores_{var_set}', sep = '\t', index = False, header = False)
    
    
    if get_tracks:
        
        if var_set == 0:
            variant_tracks_all = variant_tracks.copy()
        else:
            variant_tracks_all.update(variant_tracks)
        
        
        

    if get_maps:
        
        if var_set == 0:
            variant_maps_all = variant_maps.copy()
        else:
            variant_maps_all.update(variant_maps)
        
        
        
    var_set_list.append(var_set)
    var_set += 1
    
    
    



# Write standard output with error messages and warnings to log file
sys.stdout = std_output
log_file.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Save results
    
    
# Write sequences to fasta file
if get_seq:
    sequences_fasta = open(f'{out_file}_sequences.fa','w')
    for seq_name, sequence in sequences_all.items():
        seq_name_line = ">" + seq_name + "\n"
        sequences_fasta.write(seq_name_line)
        sequence_line = sequence + "\n"
        sequences_fasta.write(sequence_line)
    sequences_fasta.close()
    
    
# Write disruption tracks and/or predictions to array
if get_tracks:
    np.save(f'{out_file}_tracks.npy', variant_tracks_all) 
if get_maps:
    np.save(f'{out_file}_maps.npy', variant_maps_all) 


# Combine subset files into one
os.system(f'for file in {out_file}_scores_*; do cat "$file" >> {out_file}_scores && rm "$file"; done')
os.system(f'for file in {out_file}_filtered_out_*; do cat "$file" >> {out_file}_filtered_out && rm "$file"; done')



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




