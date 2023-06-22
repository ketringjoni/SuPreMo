#!/usr/bin/env python
# coding: utf-8


# # # # # # # # # # # # # # # # # # 
# # # # # Import packages # # # # #

import pandas as pd
import numpy as np

import os
import io

import math
import pysam
from scipy.stats import spearmanr
from collections import Counter

from Bio.Seq import Seq

import cooltools
from cooltools.lib.numutils import observed_over_expected
from cooltools.lib.numutils import adaptive_coarsegrain
from cooltools.lib.numutils import interpolate_bad_singletons
from cooltools.lib.numutils import interp_nan, set_diag
from cooltools.lib.plotting import *

from pathlib import Path
import gzip


svlen_limit = None
var_set_size = None
nt = ['A', 'T', 'C', 'G']



# # # # # # # # # # # # # # # # # # 
# # # # # # Load model # # # # # #

import json

from basenji import dataset
from basenji import seqnn
from basenji import dna_io
from basenji import layers

import tensorflow as tf

if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()
    
    
fasta_open = None
chrom_lengths = None
centromere_coords = None



model_file  = './Akita_model/model_best.h5'
params_file = './Akita_model/params.json'

if not Path(model_file).is_file():
    os.system('wget -P ./Akita_model/ https://storage.googleapis.com/basenji_hic/1m/models/9-14/model_best.h5')
    print('Model file downloaded as Akita_model/model_best.h5.')
if not Path(params_file).is_file():
    os.system('get -P ./Akita_model/ https://raw.githubusercontent.com/calico/basenji/master/manuscripts/akita/params.json')
    print('Model file downloaded as Akita_model/params.json.')

with open(params_file) as params_open:
    params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

params_model['augment_shift']=0

seq_length = params_model['seq_length']
target_length = params_model['target_length']

seqnn_model = seqnn.SeqNN(params_model)

hic_diags = 2
tlen = (target_length-hic_diags) * (target_length-hic_diags+1) // 2

bin_size = seq_length//target_length

seqnn_model.restore(model_file)
print('Akita successfully loaded')

hic_params = params['model']['head_hic']
cropping = hic_params[5]['cropping']
target_length_cropped = target_length - 2 * cropping

MB = None
half_patch_size = None



# # # # # # # # # # # # # # # # # # 
# # # # Get scoring code # # # # #

if not Path('scoring.py').is_file():
    os.system('wget https://raw.githubusercontent.com/pollardlab/contact_map_scoring/main/code/scoring.py')

import scoring




# # # # # # # # # # # # # # # # # # 
# # # # Reading functions # # # # #


def read_vcf(path):
    
    '''
    Read  vcf files into dataframe.
    Adapted from: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744.
    
    '''
    
    with open(path, 'r') as f:
        
        lines =[l for l in f if not l.startswith('#')]
        
    vcf_file = io.StringIO(''.join(lines))
    
    return vcf_file



def read_vcf_gz(path):
    
    '''
    Read  gzipped vcf files into dataframe.
    Adapted from: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744.
    
    '''
    
    with io.TextIOWrapper(gzip.open(path,'r')) as f:

        lines =[l for l in f if not l.startswith('#')]
        
    vcf_file = io.StringIO(''.join(lines))
    
    return vcf_file



def read_input(in_file, var_set):

    '''
    Read and reformat variant dataset. Accepted formats are .vcf .vcf.gz from 4.1 version, 
    .bed file with the following columns: [CHROM, POS, REF, ALT, END, SVTYPE, SVLEN], 
    and .tsv from ANNOVAR annotSV.
    
    '''
    
    
    if 'vcf' in in_file:
        
        # For gzipped files
        if in_file.endswith('.gz'):
            vcf_file = read_vcf_gz(in_file)
        else:
            vcf_file = read_vcf(in_file)
            
            
        variants = pd.read_csv(
                vcf_file,
                skiprows = var_set*var_set_size, nrows = var_set_size,
                sep='\t',
                names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO2', 'INFO3', 'INFO4']
            )
            
        # Read SVs
        if any(['SVTYPE' in x for x in variants.INFO]):
            
            if any(['END' in x and 'SVLEN' in x for x in variants.INFO]): # BNDs don't have 'END'
                variants['END'] = variants.INFO.str.split('END=').str[1].str.split(';').str[0] # this SVLEN (END-POS) would be 0 for SNPs
                variants.loc[~pd.isnull(variants.END), 'END'] = variants.loc[~pd.isnull(variants.END), 'END'].astype('int')
                variants['SVLEN'] = variants.INFO.str.split('SVLEN=').str[1].str.split(';').str[0]
            else:
                variants['END'] = [np.nan]*len(variants)
                variants['SVLEN'] = [np.nan]*len(variants)
            variants['SVTYPE'] = variants.INFO.str.split('SVTYPE=').str[1].str.split(';').str[0]
            
            variants = variants[['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN']]

        # Read simple variants 
        else:
            
            variants = variants[['CHROM', 'POS', 'REF', 'ALT']]       
            
 
        
    elif 'bed' in in_file:
        
        colnames = ['CHROM', 'POS', 'REF', 'ALT', 'END', 'SVTYPE', 'SVLEN']
        ncols = len(pd.read_csv(in_file, sep = '\t', nrows = 0, low_memory=False).columns)

        variants = pd.read_csv(in_file, sep = '\t', names = colnames[:ncols], low_memory=False,
                               skiprows = var_set*var_set_size, nrows = var_set_size)
        
        
    elif 'tsv' in in_file:
        
        with (gzip.open if in_file.endswith(".gz") else open)(in_file, "rt", encoding="utf-8") as variants:
            variants = (pd.read_csv(in_file, sep = '\t', low_memory=False,
                                   skiprows = var_set*var_set_size, nrows = var_set_size)
                        .rename(columns = {'SV_chrom':'CHROM', 
                                           'SV_start':'POS',
                                           'SV_end':'END', 
                                           'SV_type':'SVTYPE',
                                           'SV_length':'SVLEN'})
                       [['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN']])
            variants['CHROM'] = ['chr' + str(x) for x in variants['CHROM']]
            variants.loc[~pd.isnull(variants.END), 'END'] = variants.loc[~pd.isnull(variants.END), 'END'].astype('int')

            
    else:
        raise ValueError('Input file type not accepted. Make sure it has the right extension.')
        
        
    variants.reset_index(inplace = True, drop = True)
    
    
    return variants




# # # # # # # # # # # # # # # # # # 
# # Upstream helper functions # # #



def get_variant_position(CHR, POS, var_len, half_left, half_right):
    
    '''
    Annotate variants based on position relative to chromosome arm ends (within half the size of the prediction window of end).
    One of 5 annotations: 
    1. chrom_mid: not near chromosome arm ends
    2. chrom_start: near chromosome start
    3. chrom_centro_left: near the left end of centromere (left arm)
    4. chrom_centro_right: near the right end of centromere (right arm)
    5. chrom_end: near chromosome end.
    5. centromere: overlaps centromere. These variants will not be processed
    
    '''
    

    # Define variant position with respect to chromosome start and end

    # Get last coordinate of chromosome
    chrom_max = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max']) 
    
    # Get centromere coordinate
    centro_start = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start'])
    centro_stop = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop'])

    # If variant too close to the beginning of chromosome
    if POS - half_left <= 0: 
        var_position = "chrom_start"
        
    # If variant in centromere
    elif POS >= centro_start and POS <= centro_stop:
        var_position = "centromere"
        
    # If variant too close to left end of centromere
    elif POS + var_len - 1 + half_right > centro_start and POS - half_left < centro_stop: 
        
        if POS > centro_stop:
            var_position = "chrom_centro_right"

        elif POS < centro_start:
            var_position = "chrom_centro_left"

    # If variant too close to the end of chromosome
    elif POS + var_len - 1 + half_right > chrom_max: 
        var_position = "chrom_end"
  
    else:
        var_position = "chrom_mid"
           
    return var_position



def get_variant_type(REF, ALT):
    
    '''
    Get the type of variant from a REF and ALT allele of sequences.
    One of 4 variant types:
    1. Deletion
    2. Insertion
    3. SNP: single nucleotide variant
    4. MNP: multiple nucleotide variant
    
    '''

    # Annotate variant as one of the 6 categories below based on REF and ALT allele
    
    if len(REF) > len(ALT):
        variant_type = "Deletion"
    elif len(REF) < len(ALT):
        variant_type = "Insertion"
        
    elif len(REF) == 1 and len(ALT) ==1:
        variant_type = "SNP"
    elif len(REF) == len(ALT) and len(REF) != 1:
        variant_type = "MNP"
    
    return variant_type



def get_bin(x):
    
    '''
    Get the bin number based on the base pair number in a sequence (ex: 2500th bp is in the 2nd bin for bins of size 2048).
    Note: x is the distance to the start of the sequence, not the distance to mat_start !!!
    
    '''
    
    x_bin = math.ceil(x/bin_size) - 32
    
    return x_bin




def crop_sequence(seq, length):
    
    '''
    Remove sequence from either end to make sequece size match required sequence length.
    
    '''
    
    to_remove = (len(seq) - length)/2

    if to_remove == 0.5:
        seq = seq[1:]
    else:
        seq = seq[math.ceil(to_remove) : -math.floor(to_remove)]
        
    return seq


def from_upper_triu(vector_repr, matrix_len, num_diags):
    
    '''
    Get a matrix from a vector representatin of the upper triangle.
    
    '''
    
    z = np.zeros((matrix_len,matrix_len))
    triu_tup = np.triu_indices(matrix_len,num_diags)
    
    z[triu_tup] = vector_repr
    
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
        
    return z + z.T





# # # # # # # # # # # # # # # # # # 
# # Generating non-BND sequences # #



class get_alleles:
    
    '''
    Get reference and alternate alleles from symbolic alleles. Only supports deletions, duplications and inversions.
    
    '''
    
    def __init__(self, CHR, POS, END):
        self.CHR = CHR
        self.POS = POS
        self.END = END

    def get_alleles_DEL(self):

        REF = fasta_open.fetch(self.CHR, self.POS - 1, self.END).upper()
        ALT = REF[0]

        return REF, ALT


    def get_alleles_DUP(self):

        # Insert duplicated sequence before POS
        ALT = fasta_open.fetch(self.CHR, self.POS - 1, self.END).upper()
        REF = ALT[0]

        return REF, ALT


    def get_alleles_INV(self):

        REF = fasta_open.fetch(self.CHR, self.POS - 1, self.END).upper()
        ALT = REF[0] + str(Seq(REF[1:]).reverse_complement())

        return REF, ALT





def adjust_seq_ends(centro_start, centro_stop, chrom_max, position, adjust, shift):
           
    '''
    Get start (adjust = 1) or end (adjust = MB) of sequence for prediction based on variant position \
    with respect to chromosome arm ends (defined in get_variant_position function).
    
    '''
    
    if position == 'chrom_start':
        seq_pos = 1 + adjust + abs(shift) # 1 is added so the position is never 0. coordinates are 1-based
        
    elif position == 'chrom_centro_right':
        seq_pos = centro_stop + adjust + abs(shift)         

    elif position == 'chrom_end':
        seq_pos = chrom_max - MB + adjust - abs(shift)       

    elif position == 'chrom_centro_left':
        seq_pos = centro_start - MB + adjust - abs(shift) 
        
        
    return seq_pos




def get_sequences(CHR, POS, REF, ALT, shift, revcomp: bool):
  
    '''
    Get reference and alternate sequence for prediction from REF and ALT alleles by incorporating ALT into the reference genome.
    Requires ALT allele to be a sequence and not a symbolic allele.
    Use positive sign for a right shift and negative for a left shift.
    revcomp: Take the reverse compliment of the resulting sequence.
    
    '''

    # Get reference sequence
    
    REF_len = len(REF)

    REF_half_left = math.ceil((MB - REF_len)/2) - shift # if the REF allele is odd, shift right
    REF_half_right = math.floor((MB - REF_len)/2) + shift

    
    # Annotate whether variant position with respect to chromosome arms ends
    if len(REF) <= len(ALT):
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right)
  
    elif len(REF) > len(ALT):       
        ALT_len = len(ALT)
        ALT_half_left = math.ceil((MB - ALT_len)/2) - shift
        ALT_half_right = math.floor((MB - ALT_len)/2) + shift   
        var_position = get_variant_position(CHR, POS, ALT_len, ALT_half_left, ALT_half_right)
    

    # Get last coordinate of chromosome
    chrom_max = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max'])
    
    # Get centromere coordinates
    centro_start = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start']) + 1
    centro_stop = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop']) + 1
    
    
    # Get start and end of reference sequence
    if var_position == "chrom_mid":
        REF_start = POS - REF_half_left
        REF_stop = REF_start + MB 
    elif var_position == "centromere":
        raise ValueError('Centromeric variant')
    else:
        REF_start = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, 0, shift)
        REF_stop = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, MB, shift)
        print("Warning: Variant not centered; too close to chromosome arm ends.")
        
        
    # Get reference sequence
    REF_seq = fasta_open.fetch(CHR, REF_start - 1, REF_stop - 1).upper()


    # Error if N composition is more than 5% of sequence
    if Counter(REF_seq)['N']/MB*100 > 5:
        raise ValueError('N composition greater than 5%')



    # Error if reference sequence does not match given REF

    if var_position == "chrom_mid":
        var_rel_pos_REF = REF_half_left

    elif var_position == "chrom_start": 
        var_rel_pos_REF = POS - abs(shift) - 1

    elif var_position == "chrom_centro_right": 
        var_rel_pos_REF = POS - centro_stop - abs(shift)

    elif var_position in ["chrom_end", "chrom_centro_left"]: 
        var_rel_pos_REF = -(REF_stop - POS)


    if REF_seq[var_rel_pos_REF : var_rel_pos_REF + REF_len] != REF:
        raise ValueError('Reference allele does not match hg38.')
            
            
            
    # Error if reference sequence is not the right length      
    if len(REF_seq) != MB:
        raise ValueError('Reference sequence generated is not the right length.')





    # For SNPs, MNPs, Insertions: 
    if len(REF) <= len(ALT):

        # Create alternate sequence: change REF sequence at position from REF to ALT

        ALT_seq = REF_seq

        ALT_seq = ALT_seq[:var_rel_pos_REF] + ALT + ALT_seq[var_rel_pos_REF + REF_len:]


        var_rel_pos_ALT = var_rel_pos_REF
        
        # Chop off ends of alternate sequence if it's longer 
        if len(ALT_seq) > len(REF_seq):
            to_remove = (len(ALT_seq) - len(REF_seq))/2

            if to_remove == 0.5:
                ALT_seq = ALT_seq[1:]
                var_rel_pos_ALT = var_rel_pos_REF - 1
                
            else:
                ALT_seq = ALT_seq[math.ceil(to_remove) : -math.floor(to_remove)]
                var_rel_pos_ALT = var_rel_pos_REF - math.ceil(to_remove)
                
            


    # For Deletions
    elif len(REF) > len(ALT):


        del_len = len(REF) - len(ALT)
        
        to_add_left = math.ceil(del_len/2)
        to_add_right = math.floor(del_len/2) 

        # Get start and end of reference sequence
        if var_position == "chrom_mid":
            ALT_start = REF_start - to_add_left
            ALT_stop = REF_stop + to_add_right

        elif var_position in ["chrom_start", "chrom_centro_right"]: 
            ALT_start = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, 0, shift)
            ALT_stop = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, MB + del_len, shift)
            
        elif var_position in ["chrom_centro_left", "chrom_end"]: 
            ALT_start = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, 0 - del_len, shift)
            ALT_stop = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, MB, shift)
        
        
        
        # Get alternate sequence
        ALT_seq = fasta_open.fetch(CHR, ALT_start - 1, ALT_stop - 1).upper()
        
        
        
        # Error if alternate sequence does not match REF at POS

        if var_position == "chrom_mid":
            var_rel_pos_ALT = REF_half_left + to_add_left

        elif var_position == "chrom_start": 
            var_rel_pos_ALT = POS - abs(shift) - 1
            
        elif var_position == "chrom_centro_right": 
            var_rel_pos_ALT = POS - centro_stop - abs(shift)
            
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            var_rel_pos_ALT = -(REF_stop - POS)
                
                
        if ALT_seq[var_rel_pos_ALT : var_rel_pos_ALT + REF_len] != REF:
            raise ValueError('Sequence for the alternate allele does not match hg38 at REF position.')


    
        # Change alternate sequence to match ALT at POS

        if var_position == "chrom_mid":
            ALT_seq = ALT_seq[:var_rel_pos_ALT - 1] + ALT + ALT_seq[var_rel_pos_ALT + REF_len - 1:] 

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:var_rel_pos_ALT + 1] + ALT + ALT_seq[var_rel_pos_ALT + 1 + REF_len:]
            
        elif var_position == "chrom_centro_right": 
            ALT_seq = ALT_seq[:var_rel_pos_ALT] + ALT + ALT_seq[var_rel_pos_ALT + REF_len:]
            
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            ALT_seq = ALT_seq[:var_rel_pos_ALT] + ALT + ALT_seq[var_rel_pos_ALT + REF_len:]

            
    if len(ALT_seq) != MB:
        raise ValueError('Alternate sequence generated is not the right length.')
         
            
    # Take reverse compliment of sequence
    if revcomp:
        REF_seq, ALT_seq = [str(Seq(x).reverse_complement()) for x in [REF_seq, ALT_seq]]

        
    return REF_seq, ALT_seq, [var_rel_pos_REF, var_rel_pos_ALT]







# # # # # # # # # # # # # # # # # # 
# # # Generating BND sequences # # #



def adjust_seq_ends_BND(CHR, position, adjust, shift):
           
    '''
    Get start (adjust = 1) or end (adjust = MB) of sequence for prediction based on variant position \
    with respect to chromosome arm ends (defined in get_variant_position function).
    Different from adjust_seq_ends because it does not require centro_start, centro_stop, and chrom_max as input.
    
    '''
    
    if position == 'chrom_start':
        seq_pos = adjust + 1 + abs(shift) # 1 is added so the position is never 0. coordinates are 1-based

    elif position == 'chrom_centro_right':
        seq_pos = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop']) + adjust + abs(shift) 

    elif position == 'chrom_end':
        seq_pos = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max']) - MB + adjust - abs(shift) 

    elif position == 'chrom_centro_left':
        seq_pos = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start']) - MB + adjust - abs(shift) 
        
    return seq_pos
    


    
  
    
    
class adjust_BND_start_end:
    
    '''
    Class of functions that get genomic coordinates for seuqences for prediciton for BNDs based on the variant position with \
    respect to chromsome arm ends of both breakends.
    
    for t[p[ and ]p]t BNDs:
    - adjust_left: left breakend is not chrom_mid, right breakend is.
    - adjust_right: right breakend is not chrom_mid, left breakend is.
    - adjust_both_start: left and right breakend are chrom_start or chrom_centro_right.
    
    for [p[t BNDs:
    - adjust_left_antisense_left: left breakend is not chrom_mid, right breakend is.
    - adjust_right_antisense_left: right breakend is not chrom_mid, left breakend is.
    - adjust_start_end_antisense_left: left breakend is chrom_start or chrom_centro_right, 
                                       right breakend is chrom_end or chrom_centro_left.
    - adjust_end_start_antisense_left: left breakend is chrom_end or chrom_centro_left, 
                                       right breakend is chrom_start or chrom_centro_right.
    
    for t]p] BNDs:
    - adjust_left_antisense_right: left breakend is not chrom_mid, right breakend is.
    - adjust_right_antisense_right: right breakend is not chrom_mid, left breakend is.
    - adjust_start_end_antisense_right: left breakend is chrom_start or chrom_centro_right, 
                                        right breakend is chrom_end or chrom_centro_left.
    - adjust_end_start_antisense_right: left breakend is chrom_end or chrom_centro_left, 
                                        right breakend is chrom_start or chrom_centro_right.
    
    
    
    '''
    
    def __init__(self, CHR_left, POS_left, left_position, 
                 CHR_right, POS_right, right_position, 
                 left_start, right_end, shift):
    
        self.CHR_left = CHR_left
        self.POS_left = POS_left
        self.left_position = left_position
        self.CHR_right = CHR_right
        self.POS_right = POS_right
        self.right_position = right_position
        self.left_start = left_start
        self.right_end = right_end
        self.shift = shift


    # General 

    def ignore(self):
        return self.left_start, self.right_end


    def raise_error(self):
        raise ValueError('Cannot generate 1Mb sequence for this chromosomal rearrangement.')



    # Sense

    def adjust_left(self):

        difference = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift) - self.left_start

        self.left_start += difference
        self.right_end += difference

        return self.left_start, self.right_end


    def adjust_right(self):

        difference = adjust_seq_ends_BND(self.CHR_right, self.right_position, MB, self.shift) - self.right_end

        self.left_start += difference
        self.right_end += difference

        return self.left_start, self.right_end


    def adjust_both_start(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, MB, self.shift)

        difference = (self.POS_left - left_start) - (self.POS_right - (right_end - MB))

        if difference <=0:
            # left side is closest to start (chromosome start or right end of centromere)
            right_end -= difference

        else:
            left_start += difference

        return left_start, right_end


    def adjust_both_end(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, MB, self.shift)

        difference = (left_start + MB - self.POS_left) - (right_end - self.POS_right)

        if difference <=0:
            # left side is closest to end (chromosome end or left end of centromere)
            right_end += difference

        else:
            left_start -= difference

        return left_start, right_end




    # Antisense left


    def adjust_left_antisense_left(self):

        difference = adjust_seq_ends_BND(self.CHR_left, self.left_position, MB, self.shift) - self.left_start

        self.left_start += difference 
        self.right_end -= difference

        return self.left_start, self.right_end



    def adjust_right_antisense_left(self):

        difference = adjust_seq_ends_BND(self.CHR_right, self.right_position, MB, self.shift) - self.right_end

        self.left_start -= difference
        self.right_end += difference 

        return self.left_start, self.right_end



    def adjust_start_end_antisense_left(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, MB, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, MB, self.shift)

        difference = (self.POS_left - (left_start - MB)) - (right_end - self.POS_right)

        if difference <=0:
            # left side is closest to edge
            right_end += difference

        else:
            left_start += difference

        return left_start, right_end



    def adjust_end_start_antisense_left(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, MB, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, MB, self.shift)

        difference = (left_start - self.POS_left) - (self.POS_right - (right_end - MB))

        if difference <=0:
            # left side is closest to start (chromosome start or right end of centromere)
            right_end -= difference

        else:
            left_start -= difference

        return left_start, right_end




    # Antisense right



    def adjust_left_antisense_right(self):

        difference = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift) - self.left_start

        self.left_start += difference
        self.right_end -= difference

        return self.left_start, self.right_end


    def adjust_right_antisense_right(self):

        difference = adjust_seq_ends_BND(self.CHR_right, self.right_position, 0, self.shift) - self.right_end

        self.left_start -= difference
        self.right_end += difference

        return self.left_start, self.right_end



    def adjust_start_end_antisense_right(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, 0, self.shift)

        difference = (self.POS_left - left_start) - (MB - (self.POS_right - right_end))

        if difference <=0:
            # left side is closest to edge
            right_end += difference

        else:
            left_start += difference

        return left_start, right_end




    def adjust_end_start_antisense_right(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, 0, self.shift)

        difference = (left_start + MB - self.POS_left) - (self.POS_right - right_end)

        if difference <=0:
            # left side is closest to edge
            right_end -= difference

        else:
            left_start -= difference

        return left_start, right_end

    
    
get_BND_ALT_sense_by_pos = {'chrom_mid_chrom_mid': ['ignore', 
                                                    'ignore', 
                                                    'ignore'], # M M
                         'chrom_mid_chrom_start': ['adjust_right', 
                                                   'adjust_right_antisense_left', 
                                                   'adjust_right_antisense_right'], # M S
                         'chrom_mid_chrom_centro_left': ['adjust_right', 
                                                         'adjust_right_antisense_left', 
                                                         'adjust_right_antisense_right'], # M S
                         'chrom_mid_chrom_centro_right': ['adjust_right', 
                                                          'adjust_right_antisense_left', 
                                                          'adjust_right_antisense_right'], # M S
                         'chrom_mid_chrom_end': ['adjust_right', 
                                                 'adjust_right_antisense_left', 
                                                 'adjust_right_antisense_right'], # M S
                         'chrom_start_chrom_mid': ['adjust_left', 
                                                   'adjust_left_antisense_left', 
                                                   'adjust_left_antisense_right'], # S M
                         'chrom_start_chrom_start': ['adjust_both_start', 
                                                     'raise_error', 
                                                     'raise_error'], # S S
                         'chrom_start_chrom_centro_left': ['raise_error',
                                                           'adjust_start_end_antisense_left',
                                                           'adjust_start_end_antisense_right'], # S E
                         'chrom_start_chrom_centro_right': ['adjust_both_start', 
                                                            'raise_error', 
                                                            'raise_error'], # S S
                         'chrom_start_chrom_end': ['raise_error',
                                                   'adjust_start_end_antisense_left',
                                                   'adjust_start_end_antisense_right'], # S E
                         'chrom_centro_left_chrom_mid': ['adjust_left', 
                                                         'adjust_left_antisense_left', 
                                                         'adjust_left_antisense_right'], # E M
                         'chrom_centro_left_chrom_start': ['raise_error',
                                                           'adjust_end_start_antisense_left',
                                                           'adjust_end_start_antisense_right'], # E S
                         'chrom_centro_left_chrom_centro_left': ['adjust_both_end', 
                                                                 'raise_error', 
                                                                 'raise_error'], # E E
                         'chrom_centro_left_chrom_centro_right': ['raise_error',
                                                                  'adjust_end_start_antisense_left',
                                                                  'adjust_end_start_antisense_right'], # E S
                         'chrom_centro_left_chrom_end': ['adjust_both_end', 
                                                         'raise_error', 
                                                         'raise_error'], # E E
                         'chrom_centro_right_chrom_mid': ['adjust_left', 
                                                          'adjust_left_antisense_left', 
                                                          'adjust_left_antisense_right'], # S M
                         'chrom_centro_right_chrom_start': ['adjust_both_start', 
                                                            'raise_error', 
                                                            'raise_error'], # S S
                         'chrom_centro_right_chrom_centro_left': ['raise_error',
                                                                  'adjust_start_end_antisense_left',
                                                                  'adjust_start_end_antisense_right'], # S E
                         'chrom_centro_right_chrom_centro_right': ['adjust_both_start', 
                                                                   'raise_error', 
                                                                   'raise_error'], # S S
                         'chrom_centro_right_chrom_end': ['raise_error',
                                                          'adjust_start_end_antisense_left',
                                                          'adjust_start_end_antisense_right'], # S E
                         'chrom_end_chrom_mid': ['adjust_left', 
                                                 'adjust_left_antisense_left', 
                                                 'adjust_left_antisense_right'], # E M
                         'chrom_end_chrom_start': ['raise_error', 
                                                   'adjust_end_start_antisense_left',
                                                   'adjust_end_start_antisense_right'], # E S
                         'chrom_end_chrom_centro_left': ['adjust_both_end', 
                                                         'raise_error', 
                                                         'raise_error'], # E E
                         'chrom_end_chrom_centro_right': ['raise_error',
                                                          'adjust_end_start_antisense_left',
                                                          'adjust_end_start_antisense_right'], # E S
                         'chrom_end_chrom_end': ['adjust_both_end', 
                                                 'raise_error', 
                                                 'raise_error']} # E E


    


def get_BND_ALT_sense(CHR_left, POS_left, left_position, CHR_right, POS_right, right_position, shift):
    
    '''
    This function applies to t[p[ and ]p]t BNDs.
    
    Get sequences corresponding to the left and right side of the alternate allele of the BND and two references.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''
    
    left_start = POS_left - half_patch_size + shift
    right_end = POS_right + half_patch_size + shift
    
    
    # Adjust start and end positions based on left and right variant positions
    adjust_BND_start_end_class = adjust_BND_start_end(CHR_left, POS_left, left_position, 
                                                      CHR_right, POS_right, right_position, 
                                                      left_start, right_end, shift)
    left_start, right_end = getattr(adjust_BND_start_end_class, 
                                    get_BND_ALT_sense_by_pos[left_position+'_'+right_position][0])()        
    
    
    # Get relative position of breakend with respect to start of map
    BND_rel_pos = POS_left - left_start
    
    if BND_rel_pos != target_length_cropped/2:
        print("Warning: Variant not centered; too close to chromosome arm ends")
    
    
    # Get sequences with adjusted sequence ends
    ALT_left = fasta_open.fetch(CHR_left, left_start, POS_left).upper()

    ALT_right = fasta_open.fetch(CHR_right, POS_right - 1, right_end - 1).upper() 

    REF_for_left = fasta_open.fetch(CHR_left, left_start, left_start + MB).upper()
    REF_for_right = fasta_open.fetch(CHR_right, right_end - 1 - MB, right_end - 1).upper() 
          
    
    return ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos
    
    


def get_BND_ALT_antisense_left(CHR_left, POS_left, left_position, CHR_right, POS_right, right_position, shift):      
    
    '''
    This function applies to [p[t BNDs.
    
    Get sequences corresponding to the left and right side of the alternate allele of the BND and two references.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''

    left_start = POS_left - (- half_patch_size + shift)
    right_end = POS_right + half_patch_size + shift

    
    # Adjust start and end positions based on left and right variant positions
    adjust_BND_start_end_class = adjust_BND_start_end(CHR_left, POS_left, left_position, 
                                                      CHR_right, POS_right, right_position, 
                                                      left_start, right_end, shift)
    left_start, right_end = getattr(adjust_BND_start_end_class, 
                                    get_BND_ALT_sense_by_pos[left_position+'_'+right_position][1])()
    
    
    # Get relative position of breakend with respect to start of map
    BND_rel_pos = left_start - POS_left
    
    if BND_rel_pos != target_length_cropped/2:
        print("Warning: Variant not centered; too close to chromosome arm ends")
      
        
    # Get sequences with adjusted sequence ends
    ALT_left_revcomp = fasta_open.fetch(CHR_left, POS_left - 1, left_start - 1).upper()
    ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())

    ALT_right = fasta_open.fetch(CHR_right, POS_right - 1, right_end - 1).upper()

    REF_for_left_revcomp = fasta_open.fetch(CHR_left, left_start - 1 - MB, left_start - 1).upper() 
    REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
    REF_for_right = fasta_open.fetch(CHR_right, right_end - 1 - MB, right_end - 1).upper() 
          
    return ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos





def get_BND_ALT_antisense_right(CHR_left, POS_left, left_position, CHR_right, POS_right, right_position, shift): 
    
    '''
    This function applies to t]p] BNDs.
    
    Get sequences corresponding to the left and right side of the alternate allele of the BND and two references.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''

    left_start = POS_left - half_patch_size + shift
    right_end = POS_right - (half_patch_size + shift)

 
    # Adjust start and end positions based on left and right variant positions
    adjust_BND_start_end_class = adjust_BND_start_end(CHR_left, POS_left, left_position, 
                                                      CHR_right, POS_right, right_position, 
                                                      left_start, right_end, shift)
    left_start, right_end = getattr(adjust_BND_start_end_class, 
                                    get_BND_ALT_sense_by_pos[left_position+'_'+right_position][2])()
    
    
    # Get relative position of breakend with respect to start of map
    BND_rel_pos = POS_left - left_start
    
    if BND_rel_pos != target_length_cropped/2:
        print("Warning: Variant not centered; too close to chromosome arm ends")
      

    # Get sequences with adjusted sequence ends
    ALT_left = fasta_open.fetch(CHR_left, left_start, POS_left).upper()

    ALT_right_revcomp = fasta_open.fetch(CHR_right, right_end, POS_right).upper()
    ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

    REF_for_left = fasta_open.fetch(CHR_left, left_start, left_start + MB).upper()
    REF_for_right_revcomp = fasta_open.fetch(CHR_right, right_end, right_end + MB).upper()
    REF_for_right = str(Seq(REF_for_right_revcomp).reverse_complement())
          
    
    return ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos
            
            


def get_sequences_BND(CHR, POS, REF, ALT, shift, revcomp):
    
    '''
    This function applies all BNDs.
    
    Get sequences corresponding to the alternate BND allele (joined breakends) and two references (one from each locus that was joined). To do so, get info for the second region (CHR2, POS2), get sequences, input the alternate allele at the breakend, and crop the ends of the sequence if necessary.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''
    
    var_position = get_variant_position(CHR, POS, 0, half_patch_size - shift, half_patch_size + shift)
    
    if var_position == 'centromere':
        raise ValueError('Centromeric variant')

    if '[' in ALT:
        
        CHR2 = ALT.split(':')[0].split('[')[1]
        POS2 = int(ALT.split('[')[1].split(':')[1])
        
        var_position2 = get_variant_position(CHR2, POS2, 0, half_patch_size - shift, half_patch_size + shift)
        
        if var_position2 == 'centromere':
            raise ValueError('Centromeric variant')

        if ALT[0] in nt:

            # t[p[
            
            ALT_t = ALT.split('[')[0]
            
            ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos = get_BND_ALT_sense(CHR, POS, var_position, 
                                                                                             CHR2, POS2, var_position2, 
                                                                                             shift)



        elif ALT[0] not in nt:

            #  [p[t
            
            ALT_t = ALT.split('[')[2]

            ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos = get_BND_ALT_antisense_left(CHR2, POS2, var_position2, 
                                                                                                      CHR, POS, var_position, 
                                                                                                      shift)


    elif ']' in ALT:
        
        CHR2 = ALT.split(':')[0].split(']')[1]
        POS2 = int(ALT.split(']')[1].split(':')[1])
        
        var_position2 = get_variant_position(CHR2, POS2, 0, half_patch_size - shift, half_patch_size + shift)
        
        if var_position2 == 'centromere':
            raise ValueError('Centromeric variant')

        if ALT[0] in nt:

            # t]p]
            
            ALT_t = ALT.split(']')[0]

            ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos = get_BND_ALT_antisense_right(CHR, POS, var_position, 
                                                                                                       CHR2, POS2, var_position2, 
                                                                                                       shift)


        elif ALT[0] not in nt:

            # ]p]t

            ALT_t = ALT.split(']')[2]

            ALT_left, ALT_right, REF_for_left, REF_for_right, BND_rel_pos = get_BND_ALT_sense(CHR2, POS2, var_position2, 
                                                                                             CHR, POS, var_position, 
                                                                                             shift)


    
    if ALT_t != REF:
        
        # Since we are already including REF, remove it from ALT_t
        if ALT_t.startswith(REF):
            ALT_t = ALT_t[1:]
        elif ALT_t.endswith(REF):
            ALT_t = ALT_t[:-1]
        else:
            raise ValueError('Unexpected format: BND ALT does not include REF')
    
        # Create ALT sequence
        ALT_seq = ALT_left + ALT_t + ALT_right
        
        # chop off the sides if longer than ~1MB
        ALT_seq = crop_sequence(ALT_seq, MB)
        
    else:
        # Create ALT sequence
        ALT_seq = ALT_left + ALT_right
        
        
    if revcomp:
        REF_for_left, REF_for_right, ALT_seq = [str(Seq(x).reverse_complement()) for x in [REF_for_left, REF_for_right, ALT_seq]]

        
        
    return REF_for_left, REF_for_right, ALT_seq, [BND_rel_pos]






# # # # # # # # # # # # # # # # # # 
# # # # get_sequences_SV # # # # #


def get_sequences_SV(CHR, POS, REF, ALT, END, SVTYPE, shift, revcomp):
    
    '''
    Get reference and alternate sequences by pulling in different functions depending on the REF and ALT in put and the SV type.
    The output is a list of sequences and a relative position for the variant.
    
    For simple variants or non-BND SVs:
    sequences is [REF_seq, ALT_seq, var_rel_pos], where var_rel_pos is [var_rel_pos_REF, var_rel_pos_ALT]
    
    For BNDs:
    sequences is [REF_seq_left, REF_seq_right, ALT_seq, var_rel_pos]
    
    REF_seq is the reference sequence (left and right for each side of the breakend if BND), ALT_seq is the alternate sequence and var_rel_pos is the position of the variant (start of non-BND or breakend of BND) relative to the start of the given sequences (ex. var_rel_pos = 1000 is the 1000th bp in the sequence).
    
    '''
    
    if all([x in nt for x in REF]) and all([x in nt for x in ALT]):
        
        sequences = get_sequences(CHR, POS, REF, ALT, shift, revcomp)        
    
    elif SVTYPE in ["DEL", "DUP", "INV"]:
        
        REF, ALT = getattr(get_alleles(CHR, POS, END), f'get_alleles_{SVTYPE}')() 
    
        if len(REF) > svlen_limit or len(ALT) > svlen_limit:
            raise ValueError(f'Variant larger than {svlen_limit} bp cutoff.')
        
        sequences = get_sequences(CHR, POS, REF, ALT, shift, revcomp)
               
    elif SVTYPE == "BND":

        sequences = get_sequences_BND(CHR, POS, REF, ALT, shift, revcomp)
        
    else:
        raise ValueError('SV type not supported.')
      
    return sequences



# # # # # # # # # # # # # # # # # # 
# # # # # Running Akita # # # # # # 

def vector_from_seq(seq):
    
    '''
    Get predicted matrix from ~1 MB input sequence using Akita. 
    pred_targets is a 1D vector that corresponds to the upper tringle of a 448x448 array with the contact frequency at each 2048 bp bin corresponding to a 917,504 bp sequence (32 bins are cropped on each end from the prediction).
    
    '''
    
    seq_1hot = dna_io.dna_1hot(seq) 

    ensemble_shifts=[0,-5,5]
    sequence = tf.keras.Input(shape=(seqnn_model.seq_length, 4), name='sequence')
    sequences = [sequence]

    if len(ensemble_shifts) > 1:
        sequences = layers.EnsembleShift(ensemble_shifts)(sequences)
    sequences

    pred_targets = seqnn_model.predict(np.expand_dims(seq_1hot,0))[0,:,0]    

    return pred_targets



# # # # # # # # # # # # # # # # # # 
# # # # Processing output # # # # #



def mat_from_vector(vector_repr, matrix_len = target_length_cropped, num_diags = 2):
    
    '''
    This applies to non-BND predcitions.
    
    Get predicted matrix from Akita predictions. 
    Output is a 448x448 array with the contact frequency at each 2048 bp bin corresponding to a 917,504 bp sequence (32 bins are cropped on each end from the prediction).
    
    '''
        
    z = np.zeros((matrix_len, matrix_len))
    triu_tup = np.triu_indices(matrix_len, num_diags)
    z[triu_tup] = vector_repr
    
    for i in range(-num_diags + 1, num_diags):
        set_diag(z, np.nan, i)
        
    return z + z.T




def get_left_BND_map(pred_matrix, BND_rel_pos_map):
    
    '''
    Take upper left quarter (or more or less if shifted) of the matrix.
    
    '''
    
    left_BND_map = pred_matrix[:int(BND_rel_pos_map),
                               :int(BND_rel_pos_map)]
    
    return left_BND_map
    
    
    
def get_right_BND_map(pred_matrix, BND_rel_pos_map):
    
    '''
    Take lower right quarter (or more or less if shifted) of the matrix.
    
    '''
    
    right_BND_map = pred_matrix[int(BND_rel_pos_map):,
                                int(BND_rel_pos_map):]
    
    return right_BND_map




def assemple_BND_mats(vector_repr_L, vector_repr_R, BND_rel_pos_map, matrix_len = target_length_cropped, num_diags = 2):
    
    '''
    This applies to BND predcitions.
    
    Get predicted matrix from Akita predictions. 
    Output is a 448x448 array with the contact frequency at each 2048 bp bin corresponding to a 917,504 bp sequence (32 bins are cropped on each end from the prediction).
    
    '''
    
    z = np.zeros((matrix_len,matrix_len))

    # make df out of matrix indices and filter to get top left and bottom right parts
    indices = np.triu_indices(matrix_len, num_diags)
    indices = pd.DataFrame(np.column_stack(indices), columns = ['rows', 'cols'])

    indices_L = tuple(indices.query('cols < @BND_rel_pos_map').T.apply(np.array, axis=1))
    indices_R = tuple(indices.query('rows >= @BND_rel_pos_map').T.apply(np.array, axis=1))
    
    z[indices_L] = vector_repr_L
    z[indices_R] = vector_repr_R
    
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
        
    return z + z.T





def mask_matrices(CHR, POS, REF, ALT, REF_pred, ALT_pred, shift):

    '''
    This applied to non-BND predicted matrices.
    
    Mask reference and alternate predicted matrices based on the type of variant.
    
    '''
    
    variant_type = get_variant_type(REF, ALT)
    
    # Get last coordinate of chromosome
    chrom_max = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max']) 

    # Get centromere coordinate
    centro_start = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start'])
    centro_stop = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop'])

    
    # Insertions: Mask REF, add nans if necessary, and mirror nans to ALT
    if variant_type in ["Insertion", "Deletion"]:
        
        if variant_type == "Deletion":
            # this works the same exact way but the commands are swapped
            REF, ALT = ALT, REF
            REF_pred, ALT_pred = ALT_pred, REF_pred

        # Get REF allele sections
        REF_len = len(REF)

        REF_half_left = math.ceil((MB - REF_len)/2) - shift # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2) + shift

        
        # Get ALT allele sections
        ALT_len = len(ALT)
        
        ALT_half_left = math.ceil((MB - ALT_len)/2) - shift
        ALT_half_right = math.floor((MB - ALT_len)/2) + shift
        
        
        # Annotate whether variant is close to beginning or end of chromosome
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right)


        # Get start and end bins of REF and ALT alleles
        if var_position == "chrom_mid":
            
            var_start = get_bin(REF_half_left - 1)
            var_end = get_bin(REF_half_left - 1 + REF_len)
            
            var_start_ALT = get_bin(ALT_half_left - 1)
            var_end_ALT = get_bin(ALT_half_left - 1 + ALT_len)

        elif var_position == "chrom_start": 
            
            var_start = get_bin(POS - 1 - abs(shift))
            var_end = get_bin(POS - 1 + REF_len - abs(shift))
            
            var_start_ALT = var_start
            var_end_ALT = get_bin(POS - 1 + ALT_len - abs(shift))

        elif var_position == "chrom_centro_left": 
            
            var_start = get_bin(POS - (centro_start - MB - abs(shift)) - 1)
            var_end = get_bin(POS - (centro_start - MB - abs(shift)) - 1 + REF_len)
            
            var_start_ALT = get_bin(POS - (centro_start - MB - abs(shift)) - 1 - ALT_len)
            var_end_ALT = var_end

        elif var_position == "chrom_centro_right": 
            
            var_start = get_bin(POS - centro_stop - 1 - abs(shift))
            var_end = get_bin(POS - centro_stop - 1 + REF_len - abs(shift))
            
            var_start_ALT = var_start
            var_end_ALT = get_bin(POS - centro_stop - 1 + ALT_len - abs(shift))

        elif var_position == "chrom_end": 
            
            var_start = get_bin(POS - (chrom_max - MB - abs(shift)) - 1)
            var_end = get_bin(POS - (chrom_max - MB - abs(shift)) - 1 + REF_len)
            
            var_start_ALT = get_bin(POS - (chrom_max - MB - abs(shift)) - 1 - ALT_len)
            var_end_ALT = var_end

        elif var_position == "centromere":
            raise ValueError('Centromeric variant')


        # Mask REF map: make variant bin(s) nan and add empty bins at the variant if applicable
        
        REF_pred_masked = REF_pred.copy()

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan
  
        
        # If the ALT allele falls on more bins than the REF allele, adjust ALT allele 
            # (add nan(s) to var and remove outside bin(s))
            # Otherwise don't mask
        
        if var_end_ALT - var_start_ALT > var_end - var_start:
            
        
            # Insert the rest of the nans corresponding to the ALT allele
            to_add = (var_end_ALT - var_start_ALT) - (var_end - var_start)

            for j in range(var_start, var_start + to_add): # range only includes the first variable 
                REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 0)
                REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 1)

            # Chop off the outside of the REF matrix 
            to_remove = len(REF_pred_masked) - 448

            if var_position == "chrom_mid":
                # remove less on the left bc that's where you put one less part of the variant with odd number of bp
                REF_pred_masked = REF_pred_masked[math.floor(to_remove/2) : -math.ceil(to_remove/2), 
                                                  math.floor(to_remove/2) : -math.ceil(to_remove/2)]
                
            elif var_position in ["chrom_start", "chrom_centro_right"]: 
                # Remove all from the right
                REF_pred_masked = REF_pred_masked[: -to_remove, 
                                                  : -to_remove]

            elif var_position in ["chrom_end", "chrom_centro_left"]: 
                # Remove all from the left
                REF_pred_masked = REF_pred_masked[to_remove :, 
                                                  to_remove :]

            assert len(REF_pred_masked) == 448, 'Masked reference matrix is not the right size.'
            
            
        

        # Mask ALT map: make all nans in REF_pred also nan in ALT_pred
        
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert len(ALT_pred_masked) == 448, 'Masked alternate matrix is not the right size.'
        
        if variant_type == "Deletion":
            # Swap back
            REF_pred_masked, ALT_pred_masked = ALT_pred_masked, REF_pred_masked
        
    
    # SNPs or MNPs: Mask REF and mirror nans to ALT
    elif variant_type in ['SNP', 'MNP']:
        
        # Get REF allele sections
        REF_len = len(REF)

        REF_half_left = math.ceil((MB - REF_len)/2)  - shift # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2) + shift


        # Annotate whether variant is close to beginning or end of chromosome
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right)


        # Get start and end bins of REF and ALT alleles
        if var_position == "chrom_mid":
            
            var_start = get_bin(REF_half_left - 1)
            var_end = get_bin(REF_half_left - 1 + REF_len)

        elif var_position == "chrom_start": 
            
            var_start = get_bin(POS - 1 + abs(shift))
            var_end = get_bin(POS - 1 + REF_len + abs(shift))

        elif var_position == "chrom_centro_left": 
            
            var_start = get_bin(POS - (centro_start - MB - abs(shift)) - 1)
            var_end = get_bin(POS - (centro_start - MB - abs(shift)) - 1 + REF_len)

        elif var_position == "chrom_centro_right": 
            
            var_start = get_bin(POS - centro_stop - 1 + abs(shift))
            var_end = get_bin(POS - centro_stop - 1 + REF_len + abs(shift))

        elif var_position == "chrom_end": 
            
            var_start = get_bin(POS - (chrom_max - MB - abs(shift)) - 1)
            var_end = get_bin(POS - (chrom_max - MB - abs(shift)) - 1 + REF_len)

        elif var_position == "centromere":
            raise ValueError('Centromeric variant')
            
            
        # Mask REF map: make variant bin(s) nan 
        
        REF_pred_masked = REF_pred.copy()

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan

        
        # Mask ALT map: make all nans in REF_pred also nan in ALT_pred
        
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert len(ALT_pred_masked) == 448, 'Masked alternate matrix is not the right size.'
        
    
    
    return REF_pred_masked, ALT_pred_masked







# # # # # # # # # # # # # # # # # # 
# # Getting disruption scores # # #



class scoring_map_methods:
    
    '''
    Define the function necessary for calculating disruption scores based on map-motivated scoring methods. 
    Source: https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1.
    
    '''
    
    def __init__(self, map_a, map_b):
        self.map_a = map_a
        self.map_b = map_b

        
    # Basic methods
    
    def mse(self): # mean squared error
        return scoring.mse(self.map_a, self.map_b)
    
    def corr(self): # spearman correlation
        return scoring.spearman(self.map_a, self.map_b)
    
    def ssi(self): # structural similarity index
        return scoring.ssim_map(self.map_a, self.map_b)
    
    def scc(self): # stratum adjusted correlation
        return scoring.scc(self.map_a, self.map_b)
    
    
    # Map-motivated methods
    
    def ins(self): # insulation_track
        return scoring.vectorMethodToScalar(scoring.insulation_track, self.map_a, self.map_b, finalCompMetric = 'mse')
    
    def di(self): # DI_track
        return scoring.vectorMethodToScalar(scoring.DI_track, self.map_a, self.map_b, finalCompMetric = 'mse')
    
    def dec(self): # decay_track
        return scoring.vectorMethodToScalar(scoring.decay_track, self.map_a, self.map_b, finalCompMetric = 'mse')
    
    def tri(self): # triangle_track
        return scoring.vectorMethodToScalar(scoring.triangle_track, self.map_a, self.map_b, finalCompMetric = 'mse')
    
    def pca(self): # contact_pca_track
        return scoring.vectorMethodToScalar(scoring.contact_pca_track, self.map_a, self.map_b, finalCompMetric = 'mse')
    
    
    # Disruption tracks
    
    def mse_track(self):

        # Get vector with mean of squared differences for each row 

        sub_map = self.map_a - self.map_b

        mse_track = []
        for i in range(0,len(sub_map)):

            # Get non-nan values in ith row
            non_nan_values = np.count_nonzero(~np.isnan(sub_map[i]))

            # if there are values in ith row, get mean of squared differences
            if non_nan_values > 0:
                row_sum = np.nansum([x**2 for x in sub_map[i]])
                row_avg = row_sum/non_nan_values

            # if the whole ith row is nan, MSE is nan for ith row
            else:
                row_avg = np.nan

            # save mean of squared differences for ith row
            mse_track.append(row_avg)

        mse_track = np.array(mse_track) 

        return mse_track
    
    
    def corr_track(self):

        # Get vector with correlation score for each row

        spearman_list = []
        for i in range(0,len(self.map_a)):

            non_nan_values = np.count_nonzero(~np.isnan(self.map_b[i]))

            # if there are values in ith row, get mean of squared differences
            if non_nan_values > 0:
                spearman_results = spearmanr(self.map_a[i], self.map_b[i], axis=0, nan_policy='omit')[0] # ignores nans
                spearman_list.append(spearman_results)
            else:
                spearman_list.append(np.nan)

        corr_track = np.array(spearman_list)

        return corr_track
        




def get_scores(CHR, POS, REF, ALT, sequences, SVTYPE, scores, shift, revcomp, get_tracks, get_maps): 
    
    '''
    Get disruption scores, disruption tracks, and/or predicted maps from variants and the sequences generated from them.
    
    '''
    
    
    var_rel_pos = sequences[-1]
    
    # Error if variant position is too close to end of prediction window
    if any([int(x) <= bin_size*32 or int(x) >= MB - bin_size*32 for x in var_rel_pos]):
        raise ValueError('Variant outside prediction window after cropping')
    
    
    # Make prediction
    sequences = [x for x in sequences if type(x) == str]
    matrices = [mat_from_vector(vector) for vector in [vector_from_seq(seq) for seq in sequences]]
    
    if revcomp:
        matrices = [np.flipud(np.fliplr(x)) for x in matrices]
    
                
    # Get processed matrices to score
                     
    if len(REF) > bin_size/2 or len(ALT) > bin_size/2:

        # mask matrices
        matrices = mask_matrices(CHR, POS, REF, ALT, matrices[0], matrices[1], shift)
    
    if SVTYPE == "BND":  
        
        BND_rel_pos_map = round(var_rel_pos[0]/bin_size - 32)
        
        # Get REF and ALT vectors, excluding diagonal 
        indexes_left = np.triu_indices(BND_rel_pos_map, 2)
        indexes_right = np.triu_indices(target_length_cropped - BND_rel_pos_map, 2)

        REF_L = get_left_BND_map(matrices[0], BND_rel_pos_map)[indexes_left]
        REF_R = get_right_BND_map(matrices[1], BND_rel_pos_map)[indexes_right]
        ALT_L = get_left_BND_map(matrices[2], BND_rel_pos_map)[indexes_left]
        ALT_R = get_right_BND_map(matrices[2], BND_rel_pos_map)[indexes_right]

        matrices = (assemple_BND_mats(REF_L, REF_R, BND_rel_pos_map),
                    assemple_BND_mats(ALT_L, ALT_R, BND_rel_pos_map))


    # Get disruption score and correlation for this variant
    scores_results = {}
    
    if get_maps:
        scores_results['maps'] = [matrices[0], matrices[1]]
        
        
    for score in scores:
        scores_results[score] = getattr(scoring_map_methods(matrices[0], matrices[1]), 
                                        score)()
        
        if get_tracks and score in ['corr', 'mse']:
            scores_results[f'{score}_track'] = getattr(scoring_map_methods(matrices[0], matrices[1]), 
                                                       f'{score}_track')()
            

    return scores_results









