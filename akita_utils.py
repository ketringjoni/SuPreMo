#!/usr/bin/env python
# coding: utf-8

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

half_patch_size = 2**19
MB = 2*half_patch_size
pixel_size = 2048
bins = 448
SV_cutoff = 700000
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


model_file  = './Akita_model/model_best.h5'
params_file = './Akita_model/params.json'

if not Path(model_file).is_file():
    os.system('wget -P ./Akita_model/ https://storage.googleapis.com/basenji_hic/1m/models/9-14/model_best.h5')
if not Path(params_file).is_file():
    os.system('get -P ./Akita_model/ https://raw.githubusercontent.com/calico/basenji/master/manuscripts/akita/params.json')

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


# # # # # # # # # # # # # # # # # # 
# # # # Reading functions # # # # #


def read_vcf(path):
    
    # Read vcf files and relabel their columns
    
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})




def read_vcf_gz(path):
    
    # Read gzipped vcf files and relabel their columns
    
    with io.TextIOWrapper(gzip.open(path,'r')) as f:

        lines =[l for l in f if not l.startswith('##')]
        # Identify columns name line and save it into a dict
        dinamic_header_as_key = []
        for liness in f:
            if liness.startswith("#CHROM"):
                dinamic_header_as_key.append(liness)

    return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM':'CHROM'})




def read_input(in_file, file_format, var_type):

    # Read and reformat variant dataset
    
    if file_format == 'df':
        if var_type == 'simple':
            variants = (pd.read_csv(in_file, skiprows = 1, sep = '\t')
                        .rename(columns = {'Chromosome':'CHROM', 
                                           'Start_Position':'POS',
                                           'End_Position':'END',
                                           'Reference_Allele':'REF',
                                           'Tumor_Seq_Allele2':'ALT'})
                       [['CHROM', 'POS', 'END', 'REF', 'ALT']])
            variants['SVLEN'] = (variants.END - variants.POS).astype('int') # this SVLEN (END-POS) would be 0 for SNPs

            # Might need to do this if there are homozygous variants where Tumor_Seq_Allele1 != Reference_Allele
            # would need to only do this for those variants and concat with the rest
            # variants = pd.melt(variants, 
            #                    id_vars = ['CHROM', 'POS', 'REF'], 
            #                    value_vars = ['ALT1', 'ALT2'], 
            #                    var_name = 'allele', 
            #                    value_name='ALT')

            for i in variants[variants.ALT == '-'].index:
                variants.loc[i,'SVTYPE'] = 'DEL'
            for i in variants[variants.REF == '-'].index:
                variants.loc[i,'SVTYPE'] = 'INS'

            var_type = 'SV' # These are formatted like SVs

        elif var_type == 'SV':
            variants = (pd.read_csv(in_file, sep = '\t', low_memory=False)
                        .rename(columns = {'SV_chrom':'CHROM', 
                                           'SV_start':'POS',
                                           'SV_end':'END', 
                                           'SV_type':'SVTYPE',
                                           'SV_length':'SVLEN'})
                       [['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN']])
            variants['CHROM'] = ['chr' + str(x) for x in variants['CHROM']]
            variants.loc[~pd.isnull(variants.END), 'END'] = variants.loc[~pd.isnull(variants.END), 'END'].astype('int')


    elif file_format == 'vcf':
        if in_file.endswith('.gz'):
            variants = read_vcf_gz(in_file)
        else:
            variants = read_vcf(in_file)
            
        if var_type == 'simple':
            variants = variants[['CHROM', 'POS', 'REF', 'ALT']]

        elif var_type == 'SV':
            variants['END'] = variants.INFO.str.split('END=').str[1].str.split(';').str[0] # this SVLEN (END-POS) would be 0 for SNPs
            variants.loc[~pd.isnull(variants.END), 'END'] = variants.loc[~pd.isnull(variants.END), 'END'].astype('int')
            variants['SVTYPE'] = variants.INFO.str.split('SVTYPE=').str[1].str.split(';').str[0]
            variants['SVLEN'] = variants.INFO.str.split('SVLEN=').str[1].str.split(';').str[0]
            variants = variants[['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN']]
            
    variants.reset_index(inplace = True, drop = True)
    
    return variants
     



# # # # # # # # # # # # # # # # # # 
# # Upstream helper functions # # #



def get_variant_position(CHR, POS, var_len, half_left, half_right, chrom_lengths, centromere_coords):

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
    
    # Get the bin number based on bp number in a sequence (ex: 2500th bp is in the 2nd bin)
    # x is the distance to the start of the sequence, not the distance to mat_start !!!
    
    x_bin = math.ceil(x/pixel_size) - 32
    
    return x_bin





# # # # # # # # # # # # # # # # # # 
# # # Generating sequence # # # # #



def get_sequence(CHR, POS, REF, ALT, shift, chrom_lengths, centromere_coords, fasta_open):
  
    # Get reference and alternate sequence from REF and ALT allele using reference genome 
    # use positive sign for a right shift and negative for a left shift

    # Get reference sequence
    
    REF_len = len(REF)

    REF_half_left = math.ceil((MB - REF_len)/2) - shift # if the REF allele is odd, shift right
    REF_half_right = math.floor((MB - REF_len)/2) + shift

    
    # Annotate whether variant is near end of chromosome arms
    if len(REF) <= len(ALT): # For SNPs, MNPs, Insertions
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right, chrom_lengths, centromere_coords)
  
    elif len(REF) > len(ALT): # For Deletions        
        ALT_len = len(ALT)
        ALT_half_left = math.ceil((MB - ALT_len)/2) - shift
        ALT_half_right = math.floor((MB - ALT_len)/2) + shift   
        var_position = get_variant_position(CHR, POS, ALT_len, ALT_half_left, ALT_half_right, chrom_lengths, centromere_coords)
    

    # Get last coordinate of chromosome
    chrom_max = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max']) 
    
    # Get centromere coordinate
    centro_start = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start'])
    centro_stop = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop'])
    
    
    # Get start and end of reference sequence

    if var_position == "chrom_mid":
        REF_start = POS - REF_half_left
        REF_stop = REF_start + MB 

    elif var_position == "chrom_start": 
        REF_start = 0 + abs(shift)
        REF_stop = MB + abs(shift)
        print("Warning: Variant not centered; too close to start of chromosome.")

    elif var_position == "chrom_centro_left": 
        REF_start = centro_start - MB - abs(shift)
        REF_stop = centro_start - abs(shift)
        print("Warning: Variant not centered; too close to left end of centromere.")
        
    elif var_position == "chrom_centro_right": 
        REF_start = centro_stop + abs(shift)
        REF_stop = centro_stop + MB + abs(shift)
        print("Warning: Variant not centered; too close to right end of centromere.")

    elif var_position == "chrom_end": 
        REF_start = chrom_max - MB - abs(shift)
        REF_stop = chrom_max - abs(shift)
        print("Warning: Variant not centered; too close to end of chromosome.")
        
    elif var_position == "centromere":
        raise ValueError('Centromeric variant')

        
        
    # Get reference sequence
    REF_seq = fasta_open.fetch(CHR, REF_start, REF_stop).upper()


    # Error if Ns are more than 5% of sequence
    if Counter(REF_seq)['N']/MB*100 > 5:
        raise ValueError('N composition greater than 5%')



    # Make sure that reference sequence matches given REF

    if var_position == "chrom_mid":
        if REF_seq[(REF_half_left - 1) : (REF_half_left - 1 + REF_len)] != REF:
            raise ValueError('Reference allele does not match hg38.')

    elif var_position == "chrom_start": 
        if REF_seq[(POS - abs(shift) - 1) : (POS - abs(shift) - 1 + REF_len)] != REF:
            raise ValueError('Reference allele does not match hg38.')

    elif var_position == "chrom_centro_right": 
        POS_adj = POS - centro_stop - abs(shift)
        if REF_seq[(POS_adj - 1) : (POS_adj - 1 + REF_len)] != REF:
            raise ValueError('Reference allele does not match hg38.')

    elif var_position in ["chrom_end", "chrom_centro_left"]: 
        if REF_seq[-(REF_stop - POS + 1) : -(REF_stop - POS + 1 - REF_len)] != REF:
            raise ValueError('Reference allele does not match hg38.')


    if len(REF_seq) != MB:
            raise ValueError('Reference sequence generated is not the right length.')





    # For SNPs, MNPs, Insertions: 
    if len(REF) <= len(ALT):

        # Create alternate sequence: change REF sequence at position from REF to ALT

        ALT_seq = REF_seq

        if var_position == "chrom_mid":
            ALT_seq = ALT_seq[:(REF_half_left - 1)] + ALT + ALT_seq[(REF_half_left - 1 + REF_len):]

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:(POS - abs(shift) - 1)] + ALT + ALT_seq[(POS - abs(shift) - 1 + REF_len):]
            
        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centro_stop - abs(shift)
            ALT_seq = ALT_seq[:(POS_adj - 1)] + ALT + ALT_seq[(POS_adj - 1 + REF_len):]

        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            ALT_seq = ALT_seq[:-(REF_stop - POS + 1)] + ALT + ALT_seq[-(REF_stop - POS + 1 - REF_len):]
            
            


        # Chop off ends of alternate sequence if it's longer 
        if len(ALT_seq) > len(REF_seq):
            to_remove = (len(ALT_seq) - len(REF_seq))/2

            if to_remove == 0.5:
                ALT_seq = ALT_seq[1:]
            else:
                ALT_seq = ALT_seq[math.ceil(to_remove) : -math.floor(to_remove)]


    # For Deletions
    elif len(REF) > len(ALT):


        del_len = len(REF) - len(ALT)
        
        to_add_left = math.ceil(del_len/2)
        to_add_right = math.floor(del_len/2) 

        # Get start and end of reference sequence
        if var_position == "chrom_mid":
            ALT_start = REF_start - to_add_left
            ALT_stop = REF_stop + to_add_right

        if var_position == "chrom_start": 
            ALT_start = 0 + abs(shift)
            ALT_stop = MB + del_len + abs(shift)
            
        if var_position == "chrom_centro_left": 
            ALT_start = centro_start - MB - del_len - abs(shift)
            ALT_stop = centro_start - abs(shift)

        if var_position == "chrom_centro_right": 
            ALT_start = centro_stop + abs(shift)
            ALT_stop = centro_stop + MB + del_len + abs(shift)

        if var_position == "chrom_end": 
            ALT_start = chrom_max - MB - del_len - abs(shift)
            ALT_stop = chrom_max - abs(shift)
            
            
        # Get alternate sequence
        ALT_seq = fasta_open.fetch(CHR, ALT_start, ALT_stop).upper()
        
        
        
        # Make sure that alternate sequence matches REF at POS

        if var_position == "chrom_mid":
            if ALT_seq[(REF_half_left - 1 + to_add_left) : (REF_half_left - 1 + to_add_left + REF_len)] != REF:
                raise ValueError('Sequence for the alternate allele does not match hg38 at REF position.')

        elif var_position == "chrom_start": 
            if ALT_seq[(POS - abs(shift) - 1) : (POS - abs(shift) - 1 + REF_len)] != REF:
                raise ValueError('Sequence for the alternate allele does not match hg38 at REF position.')
            
        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centro_stop
            if ALT_seq[(POS_adj - abs(shift) - 1) : (POS_adj - abs(shift) - 1 + REF_len)] != REF:
                raise ValueError('Sequence for the alternate allele does not match hg38 at REF position.')
            
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            if ALT_seq[-(REF_stop - POS + 1) : -(REF_stop - POS - REF_len + 1)] != REF:
                raise ValueError('Sequence for the alternate allele does not match hg38 at REF position.')


    
        # Change alternate sequence to match ALT at POS

        if var_position == "chrom_mid":
            # [:N] does not include N but [N:] includes N
            ALT_seq = ALT_seq[:(REF_half_left - 1 + to_add_left)] + ALT + ALT_seq[(REF_half_left - 1 + to_add_left + REF_len):] 

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:(POS - abs(shift) - 1)] + ALT + ALT_seq[(POS - abs(shift) - 1 + REF_len):]
            
        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centro_stop
            ALT_seq = ALT_seq[:(POS_adj - abs(shift) - 1)] + ALT + ALT_seq[(POS_adj - abs(shift) - 1 + REF_len):]
            
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            ALT_seq = ALT_seq[:-(REF_stop - POS + 1)] + ALT + ALT_seq[-(REF_stop - POS - REF_len + 1):]

            
    if len(ALT_seq) != MB:
        raise ValueError('Alternate sequence generated is not the right length.')
         
        
    return REF_seq, ALT_seq





def get_sequences_BND(CHR, POS, ALT, shift, fasta_open):

    if '[' in ALT:

        if ALT[0] in nt:

            # t[p[

            CHR2 = ALT.split(':')[0].split('[')[1]
            POS2 = int(ALT.split('[')[1].split(':')[1])
            ALT_t = ALT.split('[')[0]

            ALT_left = fasta_open.fetch(CHR, POS - half_patch_size + shift, POS).upper() # don't inlcude POS

            ALT_right = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size + shift).upper() 

            REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size + shift, POS + half_patch_size + shift).upper()
            REF_for_right = fasta_open.fetch(CHR2, POS2 - half_patch_size + shift, POS2 + half_patch_size + shift).upper() 


        elif ALT[0] not in nt:

            #  [p[t

            CHR2 = ALT.split(':')[0].split('[')[1]
            POS2 = int(ALT.split('[')[1].split(':')[1])
            ALT_t = ALT.split('[')[2]

            ALT_left_revcomp = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size - shift).upper() # don't include POS2
            ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())

            ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size + shift).upper()

            REF_for_left_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size - shift, POS2 + half_patch_size - shift).upper() 
            REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
            REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size + shift, POS + half_patch_size + shift).upper() 


    elif ']' in ALT:

        if ALT[0] in nt:

            # t]p]

            CHR2 = ALT.split(':')[0].split(']')[1]
            POS2 = int(ALT.split(']')[1].split(':')[1])
            ALT_t = ALT.split(']')[0]

            ALT_left = fasta_open.fetch(CHR, POS - half_patch_size + shift, POS).upper() # don't include POS

            ALT_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size - shift, POS2).upper() # don't include POS2
            ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

            REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size + shift, POS + half_patch_size + shift).upper()
            REF_for_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size - shift, POS2 + half_patch_size - shift).upper()
            REF_for_right = str(Seq(REF_for_right_revcomp).reverse_complement())



        elif ALT[0] not in nt:

            # ]p]t

            CHR2 = ALT.split(':')[0].split(']')[1]
            POS2 = int(ALT.split(']')[1].split(':')[1])
            ALT_t = ALT.split(']')[2]

            ALT_left = fasta_open.fetch(CHR2, POS2 - half_patch_size + shift, POS2).upper() # don't include POS2

            ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size + shift).upper()

            REF_for_left = fasta_open.fetch(CHR2, POS2 - half_patch_size + shift, POS2 + half_patch_size + shift).upper() 
            REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size + shift, POS + half_patch_size + shift).upper() 


    ALT_seq = ALT_left + ALT_t + ALT_right

    # chop off the sides if longer than ~1MB
    if len(ALT_seq) > MB:
        to_remove = (len(ALT_seq) - MB)/2

        if to_remove == 0.5:
            ALT_seq = ALT_seq[1:]
        else:
            ALT_seq = ALT_seq[math.ceil(to_remove) : -math.floor(to_remove)]

        
    return REF_for_left, REF_for_right, ALT_seq








# # # # # # # # # # # # # # # # # # 
# # # # # Running Akita # # # # # # 

def vector_from_seq(seq):
    
    # Get predicted matrix from ~1MB sequence
    
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



def from_upper_triu(vector_repr, matrix_len, num_diags):
    z = np.zeros((matrix_len,matrix_len))
    triu_tup = np.triu_indices(matrix_len,num_diags)
    z[triu_tup] = vector_repr
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
    return z + z.T



def mat_from_vector(pred_targets):
    
    mat = from_upper_triu(pred_targets,target_length_cropped,2)
    #mat = interp_all_nans(mat) 

    return mat




def get_left_BND_map(pred_matrix, shift):
    # take upper left quarter (or more or less if shifted) of the matrix
    left_BND_map = pred_matrix[:int(bins/2 - round(shift/pixel_size)),
                               :int(bins/2 - round(shift/pixel_size))]
    
    return left_BND_map
    
def get_right_BND_map(pred_matrix, shift):
    # take lower right quarter (or more or less if shifted) of the matrix
    right_BND_map = pred_matrix[int(bins/2 - round(shift/pixel_size)):,
                                int(bins/2 - round(shift/pixel_size)):]
    
    return right_BND_map





def mask_matrices(CHR, POS, REF, ALT, REF_pred, ALT_pred, shift, chrom_lengths, centromere_coords):

    # Mask reference and alternate predicted matrices, based on the type of variant, when they are centered in the sequence
    
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
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right, chrom_lengths, centromere_coords)


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
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right, chrom_lengths, centromere_coords)


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


def get_MSE_from_vector(vector1, vector2):

    # Get disruption score between wt and variant vectors
    
    # A vecotr of sum of squared differences for each row of matrix 
    sub_vec = [x - y for x,y in zip(vector1, vector2)]
    
    # Get number of bins that are not nan (not including ones that overlap deletion)
    non_nan_values = np.count_nonzero(np.invert(np.isnan(sub_vec)))
    
    MSE = np.nansum([x**2 for x in sub_vec])/non_nan_values

    return MSE




def get_scores(CHR, POS, REF, ALT, score, shift, chrom_lengths, centromere_coords, fasta_open):
    
    REF_seq, ALT_seq = get_sequence(CHR, POS, REF, ALT, shift, chrom_lengths, centromere_coords, fasta_open)

    REF_vector, ALT_vector = vector_from_seq(REF_seq), vector_from_seq(ALT_seq)

    if len(REF) > pixel_size/2 or len(ALT) > pixel_size/2:

        REF_pred, ALT_pred = mat_from_vector(REF_vector), mat_from_vector(ALT_vector)

        # mask matrices
        REF_pred, ALT_pred = mask_matrices(CHR, POS, REF, ALT, REF_pred, ALT_pred, shift, chrom_lengths, centromere_coords)

        # get masked vectors
        REF_vector = REF_pred[np.triu_indices(len(REF_pred), 2)]
        ALT_vector = ALT_pred[np.triu_indices(len(ALT_pred), 2)]


    if score in ['corr', 'both']:
        correlation, corr_pval = spearmanr(REF_vector, ALT_vector, nan_policy='omit')
        if corr_pval >= 0.05:
            correlation = 1
    if score in ['mse', 'both']:
        mse = get_MSE_from_vector(REF_vector, ALT_vector)
    
    if score == 'corr':
        scores = correlation
    elif score == 'mse':
        scores = mse
    elif score == 'both':
        scores = [mse, correlation]
   
    return scores


def get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred, shift):
    
    # Get REF and ALT vectors, excluding diagonal 
    indexes_left = np.triu_indices(bins/2 - round(shift/pixel_size), 2)
    indexes_right = np.triu_indices(bins/2 + round(shift/pixel_size), 2)

    REF_L = get_left_BND_map(REF_pred_L, shift)[indexes_left]
    REF_R = get_right_BND_map(REF_pred_R, shift)[indexes_right]
    ALT_L = get_left_BND_map(ALT_pred, shift)[indexes_left]
    ALT_R = get_right_BND_map(ALT_pred, shift)[indexes_right]
    
    REF_vector = np.append(REF_L, REF_R)
    ALT_vector = np.append(ALT_L, ALT_R)
    
    # Get disruption score 
    mse = get_MSE_from_vector(REF_vector, ALT_vector)
    
    # Get spearman correlation
    correlation, corr_pval = spearmanr(REF_vector, ALT_vector)
    
    if corr_pval >= 0.05:
        correlation = 0
    
    return mse, correlation






def get_scores_SV(CHR, POS, ALT, END, SVTYPE, score, shift, chrom_lengths, centromere_coords, fasta_open):
    
    # Get new REF and ALT alleles

    if SVTYPE in ["DEL", "DUP", "INV"]:
        
        if SVTYPE == "DEL":

            # Get REF and ALT allele sequences first
            REF = fasta_open.fetch(CHR, POS - 1, END).upper()
            ALT = REF[0]


        elif SVTYPE == "DUP":

            # Insert duplicated sequence before POS
            ALT = fasta_open.fetch(CHR, POS - 1, END).upper()
            REF = ALT[0]


        elif SVTYPE == "INV":

            REF = fasta_open.fetch(CHR, POS - 1, END).upper()
            ALT = REF[0] + str(Seq(REF[1:]).reverse_complement())

        
        if len(REF) > SV_cutoff or len(ALT) > SV_cutoff:
            raise ValueError(f'Variant larger than {SV_cutoff} bp cutoff.')
            
        else:
            mse, correlation = get_scores(CHR, POS, REF, ALT, score, shift, chrom_lengths, centromere_coords, fasta_open)
     
    
        
    elif SVTYPE == "BND":

        REF_seq_L, REF_seq_R, ALT_seq = get_sequences_BND(CHR, POS, ALT, shift, fasta_open)


        REF_pred_L, REF_pred_R, ALT_pred = [mat_from_vector(vector) for vector in \
                                            [vector_from_seq(seq) for seq in [REF_seq_L, REF_seq_R, ALT_seq]]]

        # Get disruption score and correlation for this variant
        mse, correlation = get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred, shift)

    
    else:
        raise ValueError('SV type not supported.')

    
    return mse, correlation









