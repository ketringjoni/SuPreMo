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

half_patch_size = 2**19
MB = 2*half_patch_size
pixel_size = 2048
bins = 448
nt = ['A', 'T', 'C', 'G']



# # # # # # # # # # # # # # # # # # 
# # # # # # Load model # # # # # #

import json

from basenji import dataset
from basenji import seqnn
from basenji import dna_io
from basenji import layers

import tensorflow as tf
#import tensor2tensor
if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()


params_file = './Akita_model/params.json'
model_file  = './Akita_model/model_best.h5'
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
print('successfully loaded')

hic_params = params['model']['head_hic']
cropping = hic_params[5]['cropping']
target_length_cropped = target_length - 2 * cropping


# # # # # # # # # # # # # # # # # # 
# # # Read vcf file function # # #


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



# # # # # # # # # # # # # # # # # # 
# # Upstream helper functions # # #



def get_variant_position(CHR, POS, var_len, half_left, half_right, hg38_lengths, centromere_coords):

    # Define variant position with respect to chromosome start and end

    # Get last coordinate of chromosome
    chrom_max = int(hg38_lengths[hg38_lengths.CHROM == CHR[3:]]['chrom_max']) 
    
    # Get centromere coordinate
    centro_start = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start'])
    centro_stop = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop'])

    # If variant too close to the beginning of chromosome
    if POS - half_left <= 0: 
        var_position = "chrom_start"
        #print("Warning: Variant not centered; too close to start of chromosome")
        
        
    # If variant too close to left end of centromere
    elif POS + var_len - 1 + half_right > centro_start and POS - half_left < centro_stop: 
        
        if POS > centro_stop:
            var_position = "chrom_centro_right"
            #print("Warning: Variant not centered; too close to right end of centromere")
            
        elif POS < centro_start:
            var_position = "chrom_centro_left"
            #print("Warning: Variant not centered; too close to left end of centromere")
    
    # If variant in centromere
    elif POS >= centro_start and POS <= centro_stop:
        var_position = "centromere"

    # If variant too close to the end of chromosome
    elif POS + var_len - 1 + half_right > chrom_max: 
        var_position = "chrom_end"
        #print("Warning: Variant not centered; too close to end of chromosome")
        
        
    else:
        var_position = "chrom_mid"
        
        
    return var_position



def get_variant_type(REF, ALT):

    # Annotate variant as one of the 6 categories below based on REF and ALT allele
    
    if len(REF) > len(ALT) and ALT in REF:
        variant_type = "Deletion"
    elif len(REF) < len(ALT) and REF in ALT:
        variant_type = "Insertion"
        
    elif len(REF) > len(ALT) and ~(ALT in REF):
        variant_type = "Del_sub"
    elif len(REF) < len(ALT) and ~(REF in ALT):
        variant_type = "Ins_sub"
        
    elif len(REF) == 1 and len(ALT) ==1:
        variant_type = "SNP"
    elif len(REF) == len(ALT) and len(REF) != 1:
        variant_type = "MNP"
    
    return variant_type



def get_bin(x):
    
    # Get the bin number based on bp number in a sequence (ex: 2500th bp is in the 2nd bin)
    
    x_bin = math.ceil(x/pixel_size) - 32
    
    return x_bin





# # # # # # # # # # # # # # # # # # 
# # # Generating sequence # # # # #



def get_sequence(CHR, POS, REF, ALT, hg38_lengths, centromere_coords, fasta_open):
  
    # Get reference and alternate sequence from REF and ALT allele using reference genome 
    # get_sequence will give an error if N composition > 5%
    
    # Get reference sequence
    
    REF_len = len(REF)

    REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
    REF_half_right = math.floor((MB - REF_len)/2)

    # Annotate whether variant is close to beginning or end of chromosome
    var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right, hg38_lengths, centromere_coords)
    
    # Get last coordinate of chromosome
    chrom_max = int(hg38_lengths[hg38_lengths.CHROM == CHR[3:]]['chrom_max']) 
    
    # Get centromere coordinate
    centro_start = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_start'])
    centro_stop = int(centromere_coords[centromere_coords.CHROM == CHR]['centro_stop'])
    
    
    # Get start and end of reference sequence

    if var_position == "chrom_mid":
        REF_start = POS - REF_half_left
        REF_stop = REF_start + MB 

    elif var_position == "chrom_start": 
        REF_start = 0
        REF_stop = MB
        
    elif var_position == "chrom_centro_left": 
        REF_start = centro_start - MB
        REF_stop = centro_start
        
    elif var_position == "chrom_centro_right": 
        REF_start = centro_stop
        REF_stop = centro_stop + MB

    elif var_position == "chrom_end": 
        REF_start = chrom_max - MB
        REF_stop = chrom_max
        
    elif var_position == "centromere":
        raise ValueError('Centromeric variant')



    # Get reference sequence
    REF_seq = fasta_open.fetch(CHR, REF_start, REF_stop).upper()


    # Warn if Ns are more than 5% of sequence
    if Counter(REF_seq)['N']/MB*100 > 5:
        #print("Warning: N composition greater than 5%")
        raise ValueError('N composition greater than 5%')



    # Make sure that reference sequence matches given REF

    if var_position == "chrom_mid":
        assert(REF_seq[(REF_half_left - 1) : (REF_half_left - 1 + REF_len)] == REF)

    elif var_position == "chrom_start": 
        assert(REF_seq[(POS - 1) : (POS - 1 + REF_len)] == REF) 

    elif var_position == "chrom_centro_right": 
        POS_adj = POS - centro_stop
        assert(REF_seq[(POS_adj - 1) : (POS_adj - 1 + REF_len)] == REF) 

    elif var_position in ["chrom_end", "chrom_centro_left"]: 
        assert(REF_seq[-(REF_stop - POS + 1) : -(REF_stop - POS + 1 - REF_len)] == REF)


    assert(len(REF_seq) == MB)





    # For SNPs, MNPs, Insertions, or Ins_subs: 
    if len(REF) <= len(ALT):

        # Create alternate sequence: change REF sequence at position from REF to ALT

        ALT_seq = REF_seq

        if var_position == "chrom_mid":
            ALT_seq = ALT_seq[:(REF_half_left - 1)] + ALT + ALT_seq[(REF_half_left - 1 + REF_len):]

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:(POS - 1)] + ALT + ALT_seq[(POS - 1 + REF_len):]
            
        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centro_stop
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


    # For Deletions of Del_subs
    elif len(REF) > len(ALT):


        del_len = len(REF) - len(ALT)
        
        to_add_left = math.ceil(del_len/2)
        to_add_right = math.floor(del_len/2)
        

        # Get start and end of reference sequence

        if var_position == "chrom_mid":
            ALT_start = REF_start - to_add_left
            ALT_stop = REF_stop + to_add_right

        elif var_position == "chrom_start": 
            ALT_start = 0
            ALT_stop = MB + del_len
            
        elif var_position == "chrom_centro_left": 
            ALT_start = centro_start - MB - del_len
            ALT_stop = centro_start

        elif var_position == "chrom_centro_right": 
            ALT_start = centro_stop
            ALT_stop = centro_stop + MB + del_len

        elif var_position == "chrom_end": 
            ALT_start = chrom_max - MB - del_len
            ALT_stop = chrom_max


        # Get alternate sequence
        ALT_seq = fasta_open.fetch(CHR, ALT_start, ALT_stop).upper()
        
        
        
        # Make sure that alternate sequence matches REF at POS

        if var_position == "chrom_mid":
            assert(ALT_seq[(REF_half_left - 1 + to_add_left) : (REF_half_left - 1 + to_add_left + REF_len)] == REF)

        elif var_position == "chrom_start": 
            assert(ALT_seq[(POS - 1) : (POS - 1 + REF_len)] == REF)
            
        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centro_stop
            assert(ALT_seq[(POS_adj - 1) : (POS_adj - 1 + REF_len)] == REF)
            
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            assert(ALT_seq[-(REF_stop - POS + 1) : -(REF_stop - POS - REF_len + 1)] == REF)


    
        # Change alternate sequence to match ALT at POS

        if var_position == "chrom_mid":
            # [:N] does not include N but [N:] includes N
            ALT_seq = ALT_seq[:(REF_half_left - 1 + to_add_left)] + ALT + ALT_seq[(REF_half_left - 1 + to_add_left + REF_len):] 

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:(POS - 1)] + ALT + ALT_seq[(POS - 1 + REF_len):]
            
        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centro_stop
            ALT_seq = ALT_seq[:(POS_adj - 1)] + ALT + ALT_seq[(POS_adj - 1 + REF_len):]
            
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            ALT_seq = ALT_seq[:-(REF_stop - POS + 1)] + ALT + ALT_seq[-(REF_stop - POS - REF_len + 1):]

            
    assert(len(ALT_seq) == MB)
        
        
    return REF_seq, ALT_seq





def get_sequences_BND(CHR, POS, ALT, fasta_open):

    if '[' in ALT:

        if ALT[0] in nt:

            # t[p[

            CHR2 = ALT.split(':')[0].split('[')[1]
            POS2 = int(ALT.split('[')[1].split(':')[1])
            ALT_t = ALT.split('[')[0]

            ALT_left = fasta_open.fetch(CHR, POS - half_patch_size, POS).upper() # don't inlcude POS

            ALT_right = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size).upper() 

            REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper()
            REF_for_right = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper() 


        if ALT[0] not in nt:

            #  [p[t

            CHR2 = ALT.split(':')[0].split('[')[1]
            POS2 = int(ALT.split('[')[1].split(':')[1])
            ALT_t = ALT.split('[')[2]

            ALT_left_revcomp = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size).upper() # don't include POS2
            ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())

            ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size).upper()

            REF_for_left_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper() 
            REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
            REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper() 


    if ']' in ALT:

        if ALT[0] in nt:

            # t]p]

            CHR2 = ALT.split(':')[0].split(']')[1]
            POS2 = int(ALT.split(']')[1].split(':')[1])
            ALT_t = ALT.split(']')[0]

            ALT_left = fasta_open.fetch(CHR, POS - half_patch_size, POS).upper() # don't include POS

            ALT_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2).upper() # don't include POS2
            ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

            REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper()
            REF_for_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper()
            REF_for_right = str(Seq(REF_for_right_revcomp).reverse_complement())



        if ALT[0] not in nt:

            # ]p]t

            CHR2 = ALT.split(':')[0].split(']')[1]
            POS2 = int(ALT.split(']')[1].split(':')[1])
            ALT_t = ALT.split(']')[2]

            ALT_left = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2).upper() # don't include POS2

            ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size).upper()

            REF_for_left = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper() 
            REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper() 


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




def upper_left(pred_matrix):
    # take upper left quarter of the matrix
    upper_left_qt = pred_matrix[:int(bins/2),:int(bins/2)]
    
    return upper_left_qt
    
def lower_right(pred_matrix):
    # take lower right quarter of the matrix
    lower_right_qt = pred_matrix[int(bins/2):,int(bins/2):]
    
    return lower_right_qt





def mask_matrices(REF, ALT, REF_pred, ALT_pred):

    # Mask reference and alternate predicted matrices, based on the type of variant, when they are centered in the sequence
    
    variant_type = get_variant_type(REF, ALT)
    
    # Insertions: Mask REF, add nans if necessary, and mirror nans to ALT
    if variant_type in ["Insertion", "Ins_sub"]:

    # start with just the middle of the chromosome

        # Adjust reference sequence

        REF_len = len(REF)

        REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2)

        # change REF allele to nans
        var_start = get_bin(REF_half_left - 1)
        var_end = get_bin(REF_half_left - 1 + REF_len)

        REF_pred_masked = REF_pred

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan


        # Get start and end of ALT allele
        ALT_len = len(ALT)
        
        ALT_half_left = math.ceil((MB - ALT_len)/2) 
        ALT_half_right = math.floor((MB - ALT_len)/2)

        var_start_ALT = get_bin(ALT_half_left - 1)
        var_end_ALT = get_bin(ALT_half_right - 1 + ALT_len)
        
        
        
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
            to_remove = (len(REF_pred_masked) - 448)/2

            REF_pred_masked = REF_pred_masked[math.floor(to_remove) : -math.ceil(to_remove), math.floor(to_remove) : -math.ceil(to_remove)]
            # remove less on the left bc that's where you put one less part of the variant with odd number of bp

            assert(len(REF_pred_masked) == 448)
            
            
        

        # Adjust alternate sequence

        # make all nans in REF_pred also nan in ALT_pred

        # change all nan
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert(len(ALT_pred_masked) == 448)
        
    
    
    # Deletions: Mask ALT, add nans if necessary, and mirror nans to REF
    elif variant_type in ["Deletion", "Del_sub"]:
        
        ALT_len = len(ALT)

        ALT_half_left = math.ceil((MB - ALT_len)/2) # if the ALT allele is odd, shift right
        ALT_half_right = math.floor((MB - ALT_len)/2)

        # change ALT allele to nans
        var_start = get_bin(ALT_half_left - 1)
        var_end = get_bin(ALT_half_left - 1 + ALT_len)

        ALT_pred_masked = ALT_pred

        ALT_pred_masked[var_start:var_end + 1, :] = np.nan
        ALT_pred_masked[:, var_start:var_end + 1] = np.nan


        
        # Get start and end of ALT allele
        REF_len = len(REF)
        
        REF_half_left = math.ceil((MB - REF_len)/2) 
        REF_half_right = math.floor((MB - REF_len)/2)

        var_start_REF = get_bin(REF_half_left - 1)
        var_end_REF = get_bin(REF_half_right - 1 + REF_len)
        
        
        
        # If the REF allele falls on more bins than the ALT allele, adjust REF allele 
            # (add nan(s) to var and remove outside bin(s))
            # Otherwise don't mask
        
        if var_end_REF - var_start_REF > var_end - var_start:
            
        
        
            # Insert the rest of the nans corresponding to the REF allele
            to_add = (var_end_REF - var_start_REF) - (var_end - var_start)

            for j in range(var_start, var_start + to_add): # range only includes the first variable 
                ALT_pred_masked = np.insert(ALT_pred_masked, j, np.nan, axis = 0)
                ALT_pred_masked = np.insert(ALT_pred_masked, j, np.nan, axis = 1)


            # Chop off the outside of the ALT matrix 
            to_remove = (len(ALT_pred_masked) - 448)/2

            ALT_pred_masked = ALT_pred_masked[math.floor(to_remove) : -math.ceil(to_remove), math.floor(to_remove) : -math.ceil(to_remove)]
            # remove less on the left bc that's where you put one less part of the variant with odd number of bp


            assert(len(ALT_pred_masked) == 448)


        # Adjust Reference sequence

        # make all nans in ALT_pred also nan in REF_pred

        # change all nan
        ALT_pred_novalues = ALT_pred_masked.copy()

        ALT_pred_novalues[np.invert(np.isnan(ALT_pred_novalues))] = 0

        REF_pred_masked = REF_pred + ALT_pred_novalues

        assert(len(REF_pred_masked) == 448)

    
    # SNPs or MNPs: Mask REF and mirror nans to ALT
    elif variant_type in ['SNP', 'MNP']:
        
        # start with just the middle of the chromosome

        # Adjust reference sequence

        REF_len = len(REF)

        REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2)

        # change REF allele to nans
        var_start = get_bin(REF_half_left - 1)
        var_end = get_bin(REF_half_left - 1 + REF_len)

        REF_pred_masked = REF_pred

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan

        
        # Adjust alternate sequence

        # make all nans in REF_pred also nan in ALT_pred

        # change all nan
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert(len(ALT_pred_masked) == 448)
        
    
    
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




def get_scores(CHR, POS, REF, ALT, hg38_lengths, centromere_coords, fasta_open):
    
    REF_seq, ALT_seq = get_sequence(CHR, POS, REF, ALT, hg38_lengths, centromere_coords, fasta_open)

    REF_vector, ALT_vector = vector_from_seq(REF_seq), vector_from_seq(ALT_seq)

    if len(REF) > pixel_size/2 or len(ALT) > pixel_size/2:

        REF_pred, ALT_pred = mat_from_vector(REF_vector), mat_from_vector(ALT_vector)

        # mask matrices
        REF_pred, ALT_pred = mask_matrices(REF, ALT, REF_pred, ALT_pred)

        # mask vectors
        REF_vector = REF_pred[np.triu_indices(len(REF_pred), 2)]
        ALT_vector = ALT_pred[np.triu_indices(len(ALT_pred), 2)]


    correlation, corr_pval = spearmanr(REF_vector, ALT_vector, nan_policy='omit')

    if corr_pval >= 0.05:
        correlation = 1

    disruption_score = get_MSE_from_vector(REF_vector, ALT_vector)
        

    return disruption_score, correlation


def get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred):
    
    # Get REF and ALT vectors, excluding diagonal 
    indexes = np.triu_indices(bins/2, 2)
    
    REF_UL = upper_left(REF_pred_L)[indexes]
    REF_LR = lower_right(REF_pred_R)[indexes]
    ALT_UL = upper_left(ALT_pred)[indexes]
    ALT_LR = lower_right(ALT_pred)[indexes]
    
    seq_REF = np.append(REF_UL, REF_LR)
    seq_ALT = np.append(ALT_UL, ALT_LR)
    
    # Get disruption score 
    mse = get_MSE_from_vector(seq_REF, seq_ALT)
    
    # Get spearman correlation
    correlation, corr_pval = spearmanr(seq_REF, seq_ALT)
    
    if corr_pval >= 0.05:
        correlation = 0
    
    return mse, correlation






def get_scores_SV(CHR, POS, ALT, END, SVTYPE, hg38_lengths, centromere_coords, fasta_open):
    
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

        
        if len(REF) > 2*half_patch_size or len(ALT) > 2*half_patch_size:
            # Variant larger than prediction window
            mse, correlation = np.nan, np.nan
            
        else:
            mse, correlation = get_scores(CHR, POS, REF, ALT, hg38_lengths, centromere_coords, fasta_open)
     
    
        
    elif SVTYPE == "BND":

        REF_seq_L, REF_seq_R, ALT_seq = get_sequences_BND(CHR, POS, ALT, fasta_open)


        REF_pred_L, REF_pred_R, ALT_pred = [mat_from_vector(vector) for vector in \
                                            [vector_from_seq(seq) for seq in [REF_seq_L, REF_seq_R, ALT_seq]]]

        # Get disruption score and correlation for this variant
        mse, correlation = get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred)

        
    # haven't addressed big insertions (<INS> only have beginning and end of inserted sequence - ignore)

    
    return mse, correlation









