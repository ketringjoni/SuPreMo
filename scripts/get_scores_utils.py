#!/usr/bin/env python
# coding: utf-8


'''
Functions that accompany SuPreMo get_scores for scoring variants for disruption to genome folding.

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
    os.system('wget -P ./Akita_model/ https://raw.githubusercontent.com/calico/basenji/master/manuscripts/akita/params.json')
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
# # Disruption scoring methods # # #

if not Path('scripts/scoring.py').is_file():
    os.system('wget -P ./scripts/ https://raw.githubusercontent.com/pollardlab/contact_map_scoring/main/code/scoring.py')
    os.system('wget -P ./scripts/ https://raw.githubusercontent.com/pollardlab/contact_map_scoring/main/code/hicrep.py')

import sys
sys.path.insert(0, './scripts/')
import scoring



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

            non_nan_values = min(np.count_nonzero(~np.isnan(self.map_a[i])), 
                                 np.count_nonzero(~np.isnan(self.map_b[i])))

            # if there are values in ith row, get mean of squared differences
            if non_nan_values > 3: # Spearman requires at least 3 values
                spearman_results = spearmanr(self.map_a[i], self.map_b[i], axis=0, nan_policy='omit')[0] # ignores nans
                spearman_list.append(spearman_results)
            else:
                spearman_list.append(np.nan)

        corr_track = np.array(spearman_list)

        return corr_track





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


def map_from_vector(vector_repr, matrix_len = target_length_cropped, num_diags = 2):
    
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




def assemple_BND_maps(vector_repr_L, vector_repr_R, BND_rel_pos_map, matrix_len = target_length_cropped, num_diags = 2):
    
    '''
    This applies to BND predcitions.
    
    Get predicted matrix from Akita predictions. 
    Output is a 448x448 array with the contact frequency at each 2048 bp bin corresponding to a 917,504 bp sequence (32 bins are cropped on each end from the prediction).
    
    '''
    
    # Create empty matrix of NAs to place values in
    z = np.zeros((matrix_len,matrix_len))

    # make df out of matrix indices and filter to get top left and bottom right parts
    indices = np.triu_indices(matrix_len, num_diags)
    indices = pd.DataFrame(np.column_stack(indices), columns = ['rows', 'cols'])

    indices_L = tuple(indices.query('cols < @BND_rel_pos_map').T.apply(np.array, axis=1))
    indices_R = tuple(indices.query('rows >= @BND_rel_pos_map').T.apply(np.array, axis=1))
    indices_NA = tuple(indices.query('cols >= @BND_rel_pos_map & rows < @BND_rel_pos_map').T.apply(np.array, axis=1))
    
    z[indices_L] = vector_repr_L
    z[indices_R] = vector_repr_R
    z[indices_NA] = np.nan
    
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
        
    return z + z.T




def get_bin(x):
    
    '''
    Get the bin number based on the base pair number in a sequence (ex: 2500th bp is in the 2nd bin for bins of size 2048).
    Note: x is the distance to the start of the sequence, not the distance to mat_start !!!
    
    '''
    
    x_bin = math.ceil(x/bin_size) - 32
    
    return x_bin






def mask_matrices(REF_pred, ALT_pred, SVTYPE, SVLEN, var_rel_pos):

    '''
    This applied to non-BND predicted matrices.
    
    Mask reference and alternate predicted matrices based on the type of variant.
    
    '''


    if SVTYPE in ['DEL', 'DUP']:
        
        # Insertions: Add nans to reference matrix and crop ends, then mirror nans to alternate matrix
        
        if SVTYPE == 'DEL':
            
            # Deletions: Same but swapping reference and alternate values
            REF_pred, ALT_pred = ALT_pred, REF_pred
            var_rel_pos.reverse()
            
            
        # Get variant relative positions in the unmasked predicted maps
        var_start_REF = get_bin(var_rel_pos[0])

        var_start_ALT = get_bin(var_rel_pos[1])
        var_end_ALT = get_bin(var_rel_pos[1] + SVLEN)

        

        #### Mask REF map
        
        REF_pred_masked = REF_pred.copy()
        ALT_pred_masked = ALT_pred.copy()

        
        # Insert the empty bins in the reference matrix where the variant is in the alternate matrix
        to_add = var_end_ALT - var_start_ALT

        for j in range(var_start_REF, var_start_REF + to_add): # range only includes the first variable 
            REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 0)
            REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 1)

            
        # Crop the outside of the reference matrix 
        to_remove_left = var_start_REF - var_start_ALT
        to_remove_right = len(REF_pred_masked) - target_length_cropped - to_remove_left

        if to_remove_left != 0:
            REF_pred_masked = REF_pred_masked[to_remove_left:, 
                                              to_remove_left:]
        if to_remove_right != 0:
            REF_pred_masked = REF_pred_masked[:-to_remove_right, 
                                              :-to_remove_right]



        
        if SVTYPE == 'DEL':
            
            # Deletions: Swap reference and alternate values back
            REF_pred_masked, ALT_pred_masked = ALT_pred_masked, REF_pred_masked
        
        
        
    
    # Inversions: Mask REF and mirror nans to ALT
    
    elif SVTYPE == 'INV':
        
        var_start = get_bin(var_rel_pos[0])
        var_end = get_bin(var_rel_pos[0] + SVLEN)
 

        # Mask REF map: make variant bin(s) nan 
        
        REF_pred_masked = REF_pred.copy()

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan

        
        # Mask ALT map: make all nans in REF_pred also nan in ALT_pred
        
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        
        
    if len(REF_pred_masked) != target_length_cropped:
        raise ValueError('Masked reference matrix is not the right size.')    
        
        
    if len(ALT_pred_masked) != target_length_cropped:
        raise ValueError('Masked alternate matrix is not the right size.')
    
    
    return REF_pred_masked, ALT_pred_masked 



def get_masked_BND_maps(matrices, rel_pos_map):

    # Get REF and ALT vectors, excluding diagonal 
    indexes_left = np.triu_indices(rel_pos_map, 2)
    indexes_right = np.triu_indices(target_length_cropped - rel_pos_map, 2)

    REF_L = get_left_BND_map(matrices[0], rel_pos_map)[indexes_left]
    REF_R = get_right_BND_map(matrices[1], rel_pos_map)[indexes_right]

    return (assemple_BND_maps(REF_L, REF_R, rel_pos_map),matrices[2])
    
        
        
        


        
        
# # # # # # # # # # # # # # # # # # 
# # # # # # get_scores # # # # # #


def get_scores(POS, SVTYPE, SVLEN, sequences, scores, shift, revcomp, get_tracks: bool, get_maps: bool): 
    
    '''
    Get disruption scores, disruption tracks, and/or predicted maps from variants and the sequences generated from them.
    
    '''
    
    
    var_rel_pos = sequences[-1]
    rel_pos_map = get_bin(var_rel_pos[0])
    
    
    # Error if variant position is too close to end of prediction window
    if any([int(x) <= bin_size*32 or int(x) >= seq_length - bin_size*32 for x in var_rel_pos]):
        raise ValueError('Variant outside prediction window after cropping.')
    
    
    
    # Make predictions
    
    sequences = [x for x in sequences if type(x) == str]
    matrices = [map_from_vector(vector) for vector in [vector_from_seq(seq) for seq in sequences]]
    
    if revcomp:
        matrices = [np.flipud(np.fliplr(x)) for x in matrices]
    
    
    
    # Mask matrices
    
    if SVTYPE != 'BND' and abs(int(SVLEN)) > bin_size/2:

        matrices = mask_matrices(matrices[0], matrices[1], SVTYPE, abs(int(SVLEN)), var_rel_pos)

        # If masking, the relative postion on the map depends on whether it's a duplication or deletion
        # If duplication, the relative position of variant in the ALT sequence should be used
        if SVTYPE == 'DUP':
            rel_pos_map = get_bin(var_rel_pos[1])

    
    if SVTYPE == "BND":  

        matrices = get_masked_BND_maps(matrices, rel_pos_map)


        
    # Save maps, disruption scores, and tracks
    
    scores_results = {}
    
    if get_maps:
        map_start_coord = POS - var_rel_pos[0] + 32*target_length_cropped
        scores_results['maps'] = [matrices[0], matrices[1], rel_pos_map, map_start_coord]
        
        
    for score in scores:
        scores_results[score] = getattr(scoring_map_methods(matrices[0], matrices[1]), 
                                        score)()
        
        if get_tracks and score in ['corr', 'mse']:
            scores_results[f'{score}_track'] = getattr(scoring_map_methods(matrices[0], matrices[1]), 
                                                       f'{score}_track')()
            

    return scores_results









