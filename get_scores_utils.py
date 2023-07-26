#!/usr/bin/env python
# coding: utf-8


'''
Functions for scoring variants that accompany get_scores.

'''


# # # # # # # # # # # # # # # # # # 
# # # # # Import packages # # # # #

import pandas as pd
import numpy as np

import os
import io

import math
import pysam
from scipy.stats import spearmanr

import cooltools
from cooltools.lib.numutils import observed_over_expected
from cooltools.lib.numutils import adaptive_coarsegrain
from cooltools.lib.numutils import interpolate_bad_singletons
from cooltools.lib.numutils import interp_nan, set_diag
from cooltools.lib.plotting import *

from pathlib import Path




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

half_patch_size = round(seq_length/2)



# # # # # # # # # # # # # # # # # # 
# # # # Get scoring code # # # # #

if not Path('scoring.py').is_file():
    os.system('wget https://raw.githubusercontent.com/pollardlab/contact_map_scoring/main/code/scoring.py')

import scoring





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

        REF_half_left = math.ceil((seq_length - REF_len)/2) - shift # if the REF allele is odd, shift right
        REF_half_right = math.floor((seq_length - REF_len)/2) + shift

        
        # Get ALT allele sections
        ALT_len = len(ALT)
        
        ALT_half_left = math.ceil((seq_length - ALT_len)/2) - shift
        ALT_half_right = math.floor((seq_length - ALT_len)/2) + shift
        
        
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
            
            var_start = get_bin(POS - (centro_start - seq_length - abs(shift)) - 1)
            var_end = get_bin(POS - (centro_start - seq_length - abs(shift)) - 1 + REF_len)
            
            var_start_ALT = get_bin(POS - (centro_start - seq_length - abs(shift)) - 1 - ALT_len)
            var_end_ALT = var_end

        elif var_position == "chrom_centro_right": 
            
            var_start = get_bin(POS - centro_stop - 1 - abs(shift))
            var_end = get_bin(POS - centro_stop - 1 + REF_len - abs(shift))
            
            var_start_ALT = var_start
            var_end_ALT = get_bin(POS - centro_stop - 1 + ALT_len - abs(shift))

        elif var_position == "chrom_end": 
            
            var_start = get_bin(POS - (chrom_max - seq_length - abs(shift)) - 1)
            var_end = get_bin(POS - (chrom_max - seq_length - abs(shift)) - 1 + REF_len)
            
            var_start_ALT = get_bin(POS - (chrom_max - seq_length - abs(shift)) - 1 - ALT_len)
            var_end_ALT = var_end


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

        REF_half_left = math.ceil((seq_length - REF_len)/2)  - shift # if the REF allele is odd, shift right
        REF_half_right = math.floor((seq_length - REF_len)/2) + shift


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
            
            var_start = get_bin(POS - (centro_start - seq_length - abs(shift)) - 1)
            var_end = get_bin(POS - (centro_start - seq_length - abs(shift)) - 1 + REF_len)

        elif var_position == "chrom_centro_right": 
            
            var_start = get_bin(POS - centro_stop - 1 + abs(shift))
            var_end = get_bin(POS - centro_stop - 1 + REF_len + abs(shift))

        elif var_position == "chrom_end": 
            
            var_start = get_bin(POS - (chrom_max - seq_length - abs(shift)) - 1)
            var_end = get_bin(POS - (chrom_max - seq_length - abs(shift)) - 1 + REF_len)
            
            
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
    if any([int(x) <= bin_size*32 or int(x) >= seq_length - bin_size*32 for x in var_rel_pos]):
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









