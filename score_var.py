#!/usr/bin/env python

# Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020).
# Written in Python v 3.7.11

'''
usage: Akita_variant_scoring [-h] --in IN_FILE [--fa FASTA]
                             [--chrom CHROM_LENGTHS]
                             [--centro CENTROMERE_COORDS]
                             [--scores {mse,corr,ssi,scc,ins,di,dec,tri,pca} [{mse,corr,ssi,scc,ins,di,dec,tri,pca} ...]]
                             [--shift_by SHIFT_WINDOW [SHIFT_WINDOW ...]]
                             [--file OUT_FILE] [--dir OUT_DIR]
                             [--limit SVLEN_LIMIT] [--seq_len SEQ_LEN]
                             [--revcomp] [--no_revcomp] [--augment]
                             [--get_seq] [--get_tracks] [--get_maps]
                             [--no_scores] [--rows ROWS]

Pipeline for scoring variants for disruption to genome folding using Akita
(Fudenberg et. al. 2020).

optional arguments:
  -h, --help            show this help message and exit
  --in IN_FILE          Input file with variants. Accepted formats are: vcf
                        4.1 and 4.2, tsv from VEP annotation output, and BED,
                        gzipped or not for all. If BED file, must have
                        columns: REF ALT SVTYPE(SV only) SVLEN(SV only).
                        Coordinates must be 1-based.
  --fa FASTA            hg38 reference genome fasta file.
  --chrom CHROM_LENGTHS
                        File with lengths of chromosomes in hg38. Columns:
                        chromosome (ex: 1), length; no header.
  --centro CENTROMERE_COORDS
                        Centromere coordinates for hg38. Columns: chromosome
                        (ex: chr1), start, end; no header.
  --scores {mse,corr,ssi,scc,ins,di,dec,tri,pca} [{mse,corr,ssi,scc,ins,di,dec,tri,pca} ...]
                        Method(s) used to calculate disruption scores. Use
                        abbreviations as follows: mse: Mean squared error
                        corr: Spearman correlation ssi: Structural similarity
                        index measure scc: Stratum adjusted correlation
                        coefficient ins: Insulation di: Directionality index
                        dec: Contact decay tri: Triangle method pca: Principal
                        component method
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
                        700000). If seq_len length specified to be less than
                        default, limit will change to 2/3s of seq_len.
  --seq_len SEQ_LEN     Length for sequences to generate. Default value is
                        based on Akita requirement. If non-default value is
                        set, get_scores must be false
  --revcomp             Make predictions for the reverse compliment of the
                        sequence.
  --no_revcomp          Make predictions without taking the the reverse
                        compliment of the sequence.
  --augment             Get the average score with 3 augmented sequences: no
                        augmentation, -1 and 1 shift, and reverse compliment.
                        This overwrites --shift, --revcomp, and --no_revcomp
                        arguments.
  --get_seq             Save sequences for the reference and alternate alleles
                        (in fa format). If get_seq is False, get_scores must
                        be True (default). Sequence name format: {var_index}_{
                        shift}_{revcomp_annot}_{seq_index}_{var_rel_pos}.
                        var_index: input row number; shift: integer that
                        window is shifted by; revcomp_annot: present only if
                        reverse compliment of sequence was taken; seq_index:
                        index for sequences generated for that variant: 0-1
                        for non-BND reference and alternate sequences and 0-2
                        for BND left and right reference sequence and
                        alternate sequence; var_rel_pos: relative position of
                        variant in sequence: list of two for non-BND variant
                        positions in reference and alternate sequence and an
                        integer for BND breakend position in reference and
                        alternate sequences. To read fasta file:
                        pysam.Fastafile(filename).fetch(seqname, start,
                        end).upper(). To get sequence names in fasta file:
                        pysam.Fastafile(filename).references.
  --get_tracks          Save disruption score tracks: scores for each
                        predicted bin column. Only possible for mse and corr.
                        Must get disruption scores to get disruption tracks
                        (cannot specify no_scores). Saved as a dictionary in a
                        numpy file. Dictionary item name format:
                        {var_index}_{track}_{shift}_{revcomp_annot}.
                        var_index: input row number; track: disruption score
                        track specified; shift: integer that window is shifted
                        by; revcomp_annot: present only if reverse compliment
                        of sequence was taken. There is one entry per variant.
                        Each entry/variant has 2 predicitons saved as a list
                        of 2 448x448 arrays. To read dictionary in python:
                        np.load(filename, allow_pickle="TRUE").item()
  --get_maps            Save predicted contact frequency maps. Must get
                        disruption scores to get disruption tracks (cannot
                        specify no_scores). Saved as a dictionary in a numpy
                        file. Dictionary item name format:
                        {var_index}_{shift}_{revcomp_annot}. var_index: input
                        row number; shift: integer that window is shifted by;
                        revcomp_annot: present only if reverse compliment of
                        sequence was taken. There is one entry per variant.
                        Each entry/variant has 2 predicitons saved as a list
                        of 2 448x448 arrays. To read dictionary in python:
                        np.load(filename, allow_pickle="TRUE").item()
  --no_scores           Get disruption scores. If get_scores is specified as
                        False, must specify get_seq as True.
  --rows ROWS           Number of rows (variants) to read at a time from
                        input.
                        
'''



# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Parse through arguments

import argparse

parser = argparse.ArgumentParser(
                    prog = 'Akita_variant_scoring',
                    description='Pipeline for scoring variants for disruption to genome folding using Akita (Fudenberg et. al. 2020).')

parser.add_argument('--in',
                    dest = 'in_file',
                    help = 'Input file with variants. Accepted formats are: vcf 4.1 and 4.2, tsv from VEP annotation output, and BED, gzipped or not for all. If BED file, must have columns: REF ALT SVTYPE(SV only) SVLEN(SV only). Coordinates must be 1-based.', 
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

parser.add_argument('--scores',
                    dest = 'scores', 
                    nargs = '+', 
                    help = 'Method(s) used to calculate disruption scores. Use abbreviations as follows: \
                            mse: Mean squared error \
                            corr: Spearman correlation \
                            ssi: Structural similarity index measure \
                            scc: Stratum adjusted correlation coefficient \
                            ins: Insulation \
                            di: Directionality index \
                            dec: Contact decay \
                            tri: Triangle method \
                            pca: Principal component method', 
                    type = str,
                    choices = ['mse', 'corr',  'ssi', 'scc', 'ins',  'di', 'dec', 'tri', 'pca'],
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
                    help = 'Maximum structural variant length to be scored (<= 700000). If seq_len length specified to be less than default, limit will change to 2/3s of seq_len.', 
                    type = int,
                    default = 700000,
                    required = False)

parser.add_argument('--seq_len',
                    dest = 'seq_len', 
                    help = 'Length for sequences to generate. Default value is based on Akita requirement. If non-default value is set, get_scores must be false', 
                    type = int,
                    default = 1048576,
                    required = False)

parser.add_argument('--revcomp',
                    dest = 'revcomp', 
                    help = 'Make predictions for the reverse compliment of the sequence.', 
                    action='store_true', # defaults to False
                    required = False)

parser.add_argument('--no_revcomp',
                    dest = 'no_revcomp', 
                    help = 'Make predictions without taking the the reverse compliment of the sequence.', 
                    action='store_false', # defaults to True
                    required = False)

parser.add_argument('--augment',
                    dest = 'augment', 
                    help = 'Get the average score with 3 augmented sequences: no augmentation, -1 and 1 shift, and reverse compliment. This overwrites --shift, --revcomp, and --no_revcomp arguments.', 
                    action='store_true', # defaults to False
                    required = False)

parser.add_argument('--get_seq',
                    dest = 'get_seq', 
                    help = 'Save sequences for the reference and alternate alleles (in fa format). If get_seq is False, get_scores must be True (default). Sequence name format: {var_index}_{shift}_{revcomp_annot}_{seq_index}_{var_rel_pos}. var_index: input row number; shift: integer that window is shifted by; revcomp_annot: present only if reverse compliment of sequence was taken; seq_index: index for sequences generated for that variant: 0-1 for non-BND reference and alternate sequences and 0-2 for BND left and right reference sequence and alternate sequence; var_rel_pos: relative position of variant in sequence: list of two for non-BND variant positions in reference and alternate sequence and an integer for BND breakend position in reference and alternate sequences. To read fasta file: pysam.Fastafile(filename).fetch(seqname, start, end).upper(). To get sequence names in fasta file: pysam.Fastafile(filename).references.', 
                    action='store_true', # defaults to False
                    required = False)

parser.add_argument('--get_tracks',
                    dest = 'get_tracks', 
                    help = 'Save disruption score tracks: scores for each predicted bin column. Only possible for mse and corr. Must get disruption scores to get disruption tracks (cannot specify no_scores). Saved as a dictionary in a numpy file. Dictionary item name format: {var_index}_{track}_{shift}_{revcomp_annot}. var_index: input row number; track: disruption score track specified; shift: integer that window is shifted by; revcomp_annot: present only if reverse compliment of sequence was taken. There is one entry per variant. Each entry/variant has 2 predicitons saved as a list of 2 448x448 arrays. To read dictionary in python: np.load(filename, allow_pickle="TRUE").item()', 
                    action='store_true', # defaults to False
                    required = False)

parser.add_argument('--get_maps',
                    dest = 'get_maps', 
                    help = 'Save predicted contact frequency maps. Must get disruption scores to get disruption tracks (cannot specify no_scores). Saved as a dictionary in a numpy file. Dictionary item name format: {var_index}_{shift}_{revcomp_annot}. var_index: input row number; shift: integer that window is shifted by; revcomp_annot: present only if reverse compliment of sequence was taken. There is one entry per variant. Each entry/variant has 2 predicitons saved as a list of 2 448x448 arrays. To read dictionary in python: np.load(filename, allow_pickle="TRUE").item()',
                    action='store_true', # defaults to False
                    required = False)

parser.add_argument('--no_scores',
                    dest = 'get_scores', 
                    help = 'Get disruption scores. If get_scores is specified as False, must specify get_seq as True.', 
                    action='store_false', # defaults to True
                    required = False)

parser.add_argument('--rows',
                    dest = 'rows', 
                    help = 'Number of rows (variants) to read at a time from input.', 
                    type = int,
                    default = 100000,
                    required = False)


args = parser.parse_args()


in_file = args.in_file
fasta_path = args.fasta
chrom_lengths_path = args.chrom_lengths
centromere_coords_path = args.centromere_coords
scores_to_use = args.scores
shift_by = args.shift_window
out_file = args.out_file
out_dir = args.out_dir
svlen_limit = args.svlen_limit
seq_len = args.seq_len
revcomp = args.revcomp
no_revcomp = args.no_revcomp
augment = args.augment
get_seq = args.get_seq
get_tracks = args.get_tracks
get_maps = args.get_maps
get_scores = args.get_scores
var_set_size = args.rows


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Adjust inputs from arguments


# Handle argument dependencies

if seq_len != 1048576:
    get_scores = False
    if svlen_limit > 0.66*seq_len:
        svlen_limit = 0.66*seq_len
        
if not revcomp and not no_revcomp:
    raise ValueError('Either revcomp and/or no_revcomp must be True.')
if not get_seq and not get_scores:
    raise ValueError('Either get_seq and/or get_scores must be True.')
    

# Adjust shift input: Remove shifts that are outside of allowed range
max_shift = 0.4*seq_len
shift_by = [x for x in shift_by if x > -max_shift and x < max_shift]


# Adjust input for taking the reverse compliment
revcomp_decision = []

if no_revcomp:
    revcomp_decision.append(False)
if revcomp:
    revcomp_decision.append(True)
revcomp_decision_i = revcomp_decision

    
# Adjust input for taking the average score from augmented sequences
if augment:
    shift_by = [-1,0,1]
    revcomp_decision = [True, False]
#     scores_to_use = [x for x in scores_to_use if x in ['mse', 'corr']]


# Create dictionaries to save sequences, maps, and disruption score tracks, if specified
if get_seq:
    sequences = {}
    
if get_maps:
    variant_maps = {}
    if get_scores == False:
        get_scores = True
        print('Must get scores to get maps: get_scores was set to True.')
    
if get_tracks:
    variant_tracks = {}
    if get_scores == False:
        get_scores = True
        print('Must get scores to get maps: get_scores was set to True.')


    




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

# Read necessary data
fasta_open = pysam.Fastafile(fasta_path)
chrom_lengths = pd.read_table(chrom_lengths_path, header = None, names = ['CHROM', 'chrom_max'])
centromere_coords = pd.read_table(centromere_coords_path, header = None, names = ['CHROM', 'centro_start', 'centro_stop'])

utils.fasta_open = fasta_open
utils.chrom_lengths = chrom_lengths
utils.centromere_coords = centromere_coords

utils.svlen_limit = svlen_limit
utils.var_set_size = var_set_size
utils.MB = seq_len
utils.half_patch_size = round(seq_len/2)
    

    
    
    
import sys
import numpy as np

# Create log file to save standard output with error messages
# std_output = sys.stdout
# log_file = open(f'{out_file}_log','w')
# sys.stdout = log_file



    
var_set = 0
var_set_list = []
while True:
    
    # Read in variants
    variants = utils.read_input(in_file, var_set)
    if len(variants) == 0:
        break
        
    # Index input based on row number and create output with same indexes
    variants['var_index'] = list(range(var_set*var_set_size, var_set*var_set_size + len(variants)))
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
    filtered_out.var_index = filtered_out.var_index.astype('int')

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

        if 'END' in variants.columns:
            END = variant.END
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

                    sequences_i = utils.get_sequences_SV(CHR, POS, REF, ALT, END, SVTYPE, shift, revcomp)


                    if get_seq:

                        # Get relative position of variant in sequence
                        var_rel_pos = str(sequences_i[-1]).replace(', ', '_')

                        for ii in range(len(sequences_i[:-1][:3])): 
                            sequences[f'{var_index}_{shift}{revcomp_annot}_{ii}_{var_rel_pos}'] = sequences_i[:-1][ii]

                    if get_scores:

                        scores = utils.get_scores(CHR, POS, REF, ALT, sequences_i, SVTYPE, 
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

        if var_set == 0:
            variant_scores.to_csv(f'{out_file}_scores_{var_set}', sep = '\t', index = False)
        if var_set > 1:
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
# sys.stdout = std_output
# log_file.close()


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
np.save(f'{out_file}_tracks.npy', variant_tracks_all) 
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



