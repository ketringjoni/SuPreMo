#!/usr/bin/python3
"""
For questions, please contact Laura Gunsalus (laura.gunsalus@gladstone.ucsf.edu)
                           or Evonne McArthur (evonne.mcarthur@gmail.com)
"""

import math
import numpy as np
import pandas as pd
from itertools import chain

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from skimage.metrics import structural_similarity as ssim
from scipy import stats
from scipy.sparse import coo_matrix
from sklearn.decomposition import PCA
from skimage.transform import resize

from hicrep import sccByDiag
from cooltools.lib import peaks

BINS = 448  # length of side of square matrix
DIAG_OFFSET = 2  # if the diagonal is offset by a number of bins:
input_map_size = 2 ** 20


# DIAG_OFFSET = 2    DIAG_OFFSET = 1     DIAG_OFFSET = 0
#  - - a b c         - a b c d           a b c d e
#  - - - d e         - - e f g           - f g h i
#  - - - - f         - - - h i           - - j k l
#  - - - - -         - - - - j           - - - l m
#  - - - - -         - - - - -           - - - - n

# This is an example of DIAG_OFFSET = 1
#  - a b b c
#  - - - d e
#  - - - - f
#  - - - - -

# -------------------- HELPER FUNCTIONS ----------------- #

# Set lower triangle to nans
# Assert no nans

def remove_missing_points_flat(flat_a, flat_b):
    """
    Remove points in vector where either flat_a OR flat_b are NaN
    Input:
        flat_a: numpy vector
        flat_b: numpy vector
    Returns:
        flat_a: numpy vector with nan positions removed
        flat_b: numpy vector with nan positions removed
    """
    mask = np.logical_or(np.isnan(flat_a), np.isnan(flat_b))
    return flat_a[~mask], flat_b[~mask]


def fill_missing_points_map(map_a, map_b, fill=0):
    """
    Fills points in matrix where either map_a OR map_b are NaN
    Input:
        map_a: numpy array
        map_b: numpy array
        fill: value to fill missing points with (default: 0)
    Returns:
        map_a: numpy array with missing positions filled with specified value
        map_b: numpy array with missing positions filled with specified value
    """
    mask = np.logical_or(np.isnan(map_a), np.isnan(map_b))
    map_a_copy, map_b_copy = np.copy(map_a), np.copy(map_b)
    map_a_copy[mask] = fill
    map_b_copy[mask] = fill
    return map_a_copy, map_b_copy


def fill_tril(contact_map, fill=np.nan):
    """ Fill the lower triangle of a matrix with a specified value
    Input:
        map: n x n numpy array
        fill: optional value for what to fill the lower triangle with (default: nan)
    Returns:
        map_filled: n x n numpy array with lower triangle filled
    """
    map_filled = contact_map.copy()
    fill_indices = np.tril_indices(map_filled.shape[0])
    map_filled[fill_indices] = fill
    return map_filled


def spearman_1D(vector_a, vector_b):
    """
    Function to calculate the spearman correlation between two 1D arrays.

    Input:
        vector_a: 1D numpy array of length(n)
        vector_b: 1D numpy array of length(n)
    Returns:
        scalar value
    """

    vector_a, vector_b = remove_missing_points_flat(vector_a, vector_b)

    spearmanr_val, pval = stats.spearmanr(vector_a, vector_b)
    return spearmanr_val


def pearson_1D(vector_a, vector_b):
    """
    Function to calculate the pearson correlation between two 1D arrays.

    Input:
        vector_a: 1D numpy array of length(n)
        vector_b: 1D numpy array of length(n)
    Returns:
        scalar value
    """

    vector_a, vector_b = remove_missing_points_flat(vector_a, vector_b)

    pearsonr_val, pval = stats.pearsonr(vector_a, vector_b)
    return pearsonr_val


def mse_1D(vector_a, vector_b):
    """
    Function to calculate the mean squared error between two 1D arrays.

    Input:
        vector_a: 1D numpy array of length(n)
        vector_b: 1D numpy array of length(n)

    Returns:
        scalar value
    """

    vector_a, vector_b = remove_missing_points_flat(vector_a, vector_b)

    mse = np.mean(np.square(vector_a - vector_b))
    return mse


# -------------------- BASIC METHODS -------------------- #

#### MSE #####
def mse(map_a, map_b):
    """
    Mean Squared Error
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Output:
        scalar: MSE between flattened map_a and map_b
    """

    flat_a = map_a.reshape(-1)
    flat_b = map_b.reshape(-1)
    mse = mse_1D(flat_a, flat_b)

    return mse


#### SPEARMAN'S RANK CORRELATION COEF ####
def spearman(map_a, map_b):
    """
    Spearman correlation between two maps.
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Returns:
        scalar: spearmanr
    """

    flat_a = map_a.reshape(-1)
    flat_b = map_b.reshape(-1)

    spearmanr = spearman_1D(flat_a, flat_b)
    return spearmanr


#### Pearson correlation #####
def pearson(map_a, map_b):
    """
    Spearman correlation between two maps.
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Returns:
        scalar: pearsonr
    """

    flat_a = map_a.reshape(-1)
    flat_b = map_b.reshape(-1)

    pearsonr = pearson_1D(flat_a, flat_b)
    return pearsonr


#### SSI #####
# SSI ranges from -1 to 1 similar to correlation, scores closest to 0 mean more disruption
def ssim_map(map_a, map_b):
    """
    SSIM between two maps.
    SSIM ranges from -1 to 1 similar to correlation,
        scores closest to 0 mean more disruption
    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Returns:
        scalar: ssim
    """

    map_a_filled, map_b_filled = fill_missing_points_map(map_a, map_b, fill=0)

    ssim_val = ssim(map_a_filled, map_b_filled)

    return ssim_val


#### SCC #####
def scc(map_a, map_b):
    """
    HiCRep stratum adjusted correlation between two maps.
    HiCRep SCC (stratum adjusted correlation coefficient)
    Input are matrices, calculated on upper triangles,nans are converted to 0

    Input:
        map_a: n x n numpy array
        map_b: n x n numpy array
    Returns:
        scalar: scc
    """
    map_a_filled, map_b_filled = fill_missing_points_map(map_a, map_b, fill=0)

    # convert the dense matrix into a sparse matrix in COOrdinate format(required for the sccByDiag function of hicrep)
    map_a_sparse = coo_matrix(map_a_filled)
    map_b_sparse = coo_matrix(map_b_filled)
    scc = sccByDiag(map_a_sparse, map_b_sparse, map_a_sparse.shape[0])
    # scc = hicrepSCC(map_a_sparse, map_b_sparse, map_a_sparse.shape[0])
    return scc


# ---------------- MAP-MOTIVATED METHODS ---------------- #

def vectorMethodToScalar(method, map_a, map_b, finalCompMetric='all', return_tracks=False):
    """
    Wrapper for the functions that output vectors or "tracks"
    and then compare the outputs of two vectors using MSE, SpearmanR, PearsonR.
    Will default with returning a dictionary with all 3 comparison measures.
    Alternatively you can specify the final comparison measure with finalCompMetric.

    Methods that can use this include: insulation_track, DI_track, calculate_decay_track,
    triangle_track, contact_pca_track,

    Input:
        method: method to convert matrix to flattened vector track (e.g. insulation_track)
        map_a: n x n numpy array
        map_b: n x n numpy array
        finalCompMetric: either 'all', 'spearmanr', 'pearsonr', or 'mse'
        return_tracks: boolean if you want to return the 1D array tracks
    Returns:
        if finalCompMetric != 'all', returns a scalar value
        if finalCompMetric == 'all', returns a dictionary of comparisons with all 3
        (spearmanr, pearsonr, and mse)
    """

    map_a_filled, map_b_filled = fill_missing_points_map(map_a, map_b, fill=np.nan)

    a_track = method(map_a_filled)
    b_track = method(map_b_filled)

    if finalCompMetric == 'all':
        output = {'spearmanr': spearman_1D(a_track, b_track),
                  'pearsonr': pearson_1D(a_track, b_track),
                  'mse': mse_1D(a_track, b_track)}
    elif finalCompMetric == 'spearmanr':
        output = spearman_1D(a_track, b_track)
    elif finalCompMetric == 'pearsonr':
        output = pearson_1D(a_track, b_track)
    elif finalCompMetric == 'mse':
        output = mse_1D(a_track, b_track)
    else:
        raise ValueError(
            f"finalCompMetric specified ({finalCompMetric}) is not understood. Use either 'all', 'spearmanr', 'pearsonr', or 'mse'")

    if return_tracks:
        return output, (a_track, b_track)

    return output


#### INSULATION TRACK ####
# Calculated on upper triangles, nans are kept but nanmean are calculated
def insulation_track(map, window_size=10, plot=False, ax=None):
    """
    get the insulation profile using a diamond-shaped window-based method, specially it scans along
    the diagonal of a matrix using a W by W diamond-shaped window, calculating the average contact
    frequency within each window. The locations at which the average contact frequency reaches a
    local minimum are identified as candidate TAD boundaries.
    Input:
        map: n x n numpy array
        window_size: size of the diamond-shaped window
        plot: True or False
        ax: if you want to plot on an already specified axis
    Returns:
        array: n x 1 insulation track of map
    """

    insulation_track = []

    # Select the diagonal
    for loc in range(0, len(map)):
        # Ignore if it's an edge pixel
        if loc <= window_size or loc >= len(map) - window_size:
            insulation_track.append(np.nan)
            continue

        # Define focal region
        focal_start = loc - window_size - 1
        focal_end = loc + window_size + 1

        window = map[focal_start:loc - 1, loc + 1:focal_end]
        window_mean = np.nanmean([np.exp(i) for i in list(chain(*window))])

        insulation_track.append(math.log2(window_mean))

    if plot:
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(insulation_track)
        ax.set_xlim(0, BINS)
        ax.set_ylabel('insulation score')

    return np.array(insulation_track)


#### CONTACT DECAY TRACK ####

def decay_track(contact_map):
    """
    Calculate contact decay track.
    Input:
        map: n x n numpy array
        plot: True or False
        ax: if you want to plot on an already specified axis
    Output:
        n x 1 contact decay track.
    """

    valid_map = np.nan_to_num(contact_map)  # replace NAs with 0s

    # get all indices
    XX, YY = np.meshgrid(np.arange(valid_map.shape[1]), np.arange(valid_map.shape[0]))

    #  get table in form value, index, index
    # (eg. 49, 0, 1 -> intensity of 49 at Ai,j, where i=0 and j=1)
    table = np.vstack((valid_map.ravel(), XX.ravel(), YY.ravel())).T

    combos = pd.DataFrame(table, columns=['val', 'x', 'y'])
    combos['distance'] = np.abs(combos['x'] - combos['y'])
    combos = combos[combos['val'] != 0]  # mask zeros

    combos['val'] = combos['val'] + abs(np.min(combos['val']))  # make positive
    combos = combos.sort_values(by='distance')  # sort values

    reduced_median = combos.groupby('distance').agg('median')

    return reduced_median['val'].values


#### TRIANGLE TRACK ####

def triangle_track(contact_map, plot=False, ax=None, plottingParams=(.25, 0.15, 30)):
    """
    Creates a track that calculates all sub-triangles of the upper triangle map
    For each sub-triangle, calculates an average contact. It
    outputs the resulting vector of means of each sub-triangles.

    For example, with a 5x5 map that looked like:
    - - a b c
    - - - d e
    - - - - f
    - - - - -
    The sub-triangles considered would be:
    a b | d e | a b c
      d |   f |   d e
        |     |     f
    Input:
        map: n x n numpy array
        plot: true or false
        ax: axis handle if specified
        plottingParams = (vmin/vmax for coloring plot, default = 0.25
                      spacing between lines, default = 0.15,
                      number of lines on plot, default = 30 [max = 445])
    Returns:
        triangle track: n x 1 numpy array
    """
    map_filled = fill_tril(contact_map.copy(), np.nan)  # set bottom triangle to nan
    map_filled_trimmed = map_filled[:-DIAG_OFFSET, DIAG_OFFSET:]  # Remove diagonal offset

    # Iterate through sub-triangles and take mean, store as vector
    length = len(map_filled_trimmed)
    triangle_track = np.array(
        [np.nanmean(map_filled_trimmed[r:, :c + 1]) for r in range(length) for c in range(length) if c > r])

    if not (plot):
        return triangle_track

    if ax is None:
        fig, ax = plt.subplots()

    plotting_indices = np.array([[r, c] for r in range(length) for c in range(length) if c > r])

    plot_df = pd.DataFrame(
        {'row': plotting_indices[:, 0], 'col': plotting_indices[:, 1], 'val': triangle_track}).set_index(
        ['row', 'col'])
    norm = plt.Normalize(-plottingParams[0],
                         plottingParams[0])  # can change this to control the dynamic range of the colors
    spacing = plottingParams[1]  # can change this to control the spacing of lines
    space = 0
    for diag in np.linspace(1, length - 1, min(plottingParams[2], length - 1)):
        plot_row_indices = np.where(np.eye(length, k=int(diag)) == 1)  # indices along each diagonal
        diag_vals = plot_df.loc[[(r, c) for r, c in zip(plot_row_indices[0], plot_row_indices[
            1])]].values.flatten()  # values corresponding to the triangle track along this diagonal

        x = np.arange(len(diag_vals)) - len(diag_vals) / 2
        y = diag_vals + space
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, cmap='RdBu_r', norm=norm)
        lc.set_array(y - space)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        plt.scatter(x, y, s=0, edgecolor='none')
        space += spacing

    return triangle_track


#### Eigenvector difference ####
def contact_pca_track(contact_map):
    """
    Return first principal component of map.

    Input:
        n x n matrix
    Output:
        n x 1 compartment track.
    """
    valid_map = np.nan_to_num(contact_map, copy=True, nan=0)  # replace NAs with 0s
    contact_model = PCA(1)
    contact_model.fit(valid_map)
    contact_pc1 = np.reshape(contact_model.fit_transform(valid_map), -1)
    return contact_pc1


#### DIRECTIONALITY INDEX TRACK ####

def downres(input_map, new_resolution=40000, input_map_size=2 ** 20):
    """
    Change the resolution of the contact map.
    input:
        input_map: n x n numpy array
        new_resolution (in bp)
        input_map_size: original map size (in bp)
    output:
        m x m numpy array, where m is the no of bp per pixel
    """
    pixel_size = round(input_map_size / new_resolution)
    resized = resize(input_map, (pixel_size, pixel_size), anti_aliasing=False)
    return resized


def DI_track(input_mat, input_map_size=2 ** 20,
             new_resolution=2000,
             window_resolution=10000,
             replace_ends=True,
             buffer=50):
    """
    Calculate directionality index track
    Source: https://zhonglab.gitbook.io/3dgenome/chapter2-computational-analysis/3.2-higer-order-data-analysis/tad-calling-algorithms

    Input:
        input_mat: n x n numpy array
        input_map_size: size of original map in bp
        new_resolution: resolution of intended map in bp
        window_resolution: resolution of sliding window in bp
        replace_ends: replaces ends of DI track with 0s
        buffer: how far to replace with 0
    Output:
        DI track
    """
    downres_map = downres(input_mat, new_resolution)
    bp_per_pixels = np.around(input_map_size / new_resolution)
    pixels_per_window = round(window_resolution / bp_per_pixels)

    summed_map = np.nansum(downres_map, axis=0)  # contact summed across one axis
    extended_map = np.concatenate([np.repeat(summed_map[0], pixels_per_window),
                                   summed_map,
                                   np.repeat(summed_map[0], pixels_per_window)])  # extend in window size each direction

    DI = []
    for i in range(pixels_per_window, summed_map.shape[0] + pixels_per_window):
        A = extended_map[i - pixels_per_window:i].sum()
        B = extended_map[i:i + pixels_per_window].sum()
        E = (A + B) / 2

        sign = (B - A) / abs(B - A)
        upstream = ((A - E) ** 2) / E
        downstream = ((B - E) ** 2) / E
        score = sign * (upstream + downstream)
        DI.append(score)
    if replace_ends:
        return np.array([0] * buffer + DI[buffer:len(DI) - buffer] + [0] * buffer)
    return np.array(DI)


# ------------ BIOLOGICALLY-INFORMED METHODS ------------ #

# ### Loop caller ####
# #input are matrices, calculated on upper triangles,nans are converted to 0 after get the exp
# values of the matrix as the the current matrices are on log scale #haven't modified the code to output the ratios,
# but if preferred after reviewing the code, could make the changes.

def findloops(matrix, p=2, width=5, ther=1.1, ther_H=1.1, ther_V=1.1):
    """
    Call loops from maps
    Input:
        matrix: n x n numpy array
        p:  the width of the interaction region surrounding the peak
        width: the size to get the donut filter
        ther: the threshold for the ratio of center windows to the donut filter and lower left filter
        ther_H: the threshold for the ratio of center windows to the horizontal filter
        ther_V: the threshold for the ratio of center windows to the vertical filter
    Return:
        loops coords
    """
    R, C = np.triu_indices(matrix.shape[0], 2)  # DIAG_OFFSET (?)
    coords = [(r, c) for r, c in zip(R, C)]
    M = np.nan_to_num(np.exp(matrix))
    Ridx = []
    Cidx = []
    w = width
    for x, y in coords:
        if (x - w < 0) or (x + w + 1 > M.shape[0]) or (y - w < 0) or (y + w + 1 > M.shape[0]):
            continue

        window = M[x - w:x + w + 1, y - w:y + w + 1]
        center = window[w, w]
        if (center < 0) or (center < np.max(window)):
            continue

        center_w = M[x - p:x + p + 1, y - p:y + p + 1]

        # lower left neighborhood
        LL_f = M[x + 1:x + w + 1, y - w:y]
        LL_c = M[x + 1:x + p + 1, y - p:y]
        # vertical filter above the pixel
        VA = M[x - w:x - p, y - 1:y + 2]
        # vertical filter below the pixel
        VB = M[x + p + 1:x + w + 1, y - 1:y + 2]
        # horizontal filter left to the pixel
        HL = M[x - 1:x + 2, y - w:y - p]
        # horizontal filter right to the pixel
        HR = M[x - 1:x + 2, y + p + 1:y + w + 1]

        D_mean = (np.sum(window) - np.sum(center_w)) / (
                window.shape[0] * window.shape[1] - center_w.shape[0] * center_w.shape[1])
        LL_mean = (np.sum(LL_f) - np.sum(LL_c)) / (LL_f.shape[0] * LL_f.shape[1] - LL_c.shape[0] * LL_c.shape[1])
        V_mean = np.mean(np.vstack((VA, VB)))
        H_mean = np.mean(np.hstack((HL, HR)))

        ZD = np.mean(center_w) / D_mean
        ZL = np.mean(center_w) / LL_mean
        ZV = np.mean(center_w) / V_mean
        ZH = np.mean(center_w) / H_mean

        if (ZD > ther) and (ZL > ther) and (ZV > ther_V) and (ZH > ther_H) and ((y - x) > 2 * width):
            Ridx.append(x)
            Cidx.append(y)
    return Ridx, Cidx


# get the same and different loops of two maps
def loops_diff(wt_matrix, del_matrix, p=2, width=5, ther=1.1, ther_H=1.1, ther_V=1.1, radius=5, detail=False):
    """
    get the same and different loops of two maps
    Input:
        wt_matrix: n x n numpy array
        del_matrix: n x n numpy array
        p:  the width of the interaction region surrounding the peak
        width: the size to get the donut filter
        ther: the threshold for the ratio of center windows to the donut filter and lower left filter
        ther_H: the threshold for the ratio of center windows to the horizontal filter
        ther_V: the threshold for the ratio of center windows to the vertical filter
        radius: the upper bound of distance of two loop points considered as same
        detail: either detailed loop information or the number of loops (all, overlap, gain, loss)
    Returns:
        if detail = True, return a dictionary of the detailed loop information of all loops in wt, del,
            overlapped loops between wt and del, loops lost in del, and loops gained in del
        if detail = False, return a dictionary of the number of all loops in wt, del, overlapped loops between
        wt and del,loops lost in del, and loops gained in del
    """
    wt_r, wt_c = findloops(wt_matrix, p=2, width=5, ther=ther, ther_H=ther_H, ther_V=ther_V)
    del_r, del_c = findloops(del_matrix, p=2, width=5, ther=ther, ther_H=ther_H, ther_V=ther_V)

    # get the overlap between loops in wt and del map
    overlap_wt = set()
    overlap_del = set()
    wt_all = set(zip(wt_r, wt_c))
    del_all = set(zip(del_r, del_c))
    for i, j in del_all:
        for k, v in wt_all:
            if abs(k - i) <= radius and abs(v - j) <= radius:
                overlap_wt.add((k, v))
                overlap_del.add((i, j))
    loss = wt_all - overlap_wt
    gain = del_all - overlap_del

    if detail:
        return {'wt': list(wt_all),
                'del': list(del_all),
                'overlap_loop_in_wt': list(overlap_wt),
                'overlap_loop_in_del': list(overlap_del),
                'loss_loop': list(loss),
                'gain_loop': list(gain)}
    else:
        return {
            'wt': len(wt_all),
            'del': len(del_all),
            'overlap': len(overlap_del),
            'overlap_ratio': len(overlap_del) / max(len(del_all), len(wt_all)),
            'gain': len(gain),
            'gain_ratio': len(gain) / max(len(del_all), len(wt_all)),
            'loss': len(loss),
            'loss_ratio': len(loss) / max(len(del_all), len(wt_all))}


# get the TAD boundaries and different between two maps
#### TADs ####
def TAD(wt_matrix, del_matrix, window_size=5, ther=0.2, radius=5, plot=False, detail=False):
    """
    identity TAD boundaries from insulation profile and get the TAD boundaries difference between two maps
    Input:
        wt_matrix: n x n numpy array
        del_matrix: n X n numpy array
        window_size: size of the diamond-shaped window
        ther: the threshold for TAD boundaries
        radius: the upper bound of distance of two TADs considered as same
        plot: True or False
        detail: print TAD boundary loci or the number of overlapped TAD boundaries
    Returns:
        if detail = True, return a dictionary of detailed TAD boundaries of wt, del, overlap, lost in del, gained in del
        if detail = False, return a dictionary of the number of TAD boundaries of wt, del, overlap,
            lost in del, gained in del
    """
    insul_track = insulation_track(wt_matrix, window_size=window_size, plot=plot)
    insul_track_del = insulation_track(del_matrix, window_size=window_size, plot=plot)

    ##  find peak prominence, calculate insulation strength
    poss, proms = peaks.find_peak_prominence(-np.array(insul_track))
    poss_del, proms_del = peaks.find_peak_prominence(-np.array(insul_track_del))
    mask = proms >= ther
    mask_del = proms_del >= ther

    overlap_wt = [x for x in poss[mask] for y in poss_del[mask_del] if abs(x - y) < radius]
    overlap_del = [y for x in poss[mask] for y in poss_del[mask_del] if abs(x - y) < radius]

    loss = list(set(poss[mask]) - set(overlap_wt))
    gain = list(set(poss_del[mask_del]) - set(overlap_del))

    if detail:
        return {'wt': list(zip(poss[mask], proms[mask])),
                'del': list(zip(poss_del[mask_del], proms_del[mask_del])),
                'overlap_TAD_boundaries_in_wt': overlap_wt,
                'overlap_TAD_boundaries_in_del': overlap_del,
                'loss_TAD': loss,
                'gain_TAD': gain}
    else:
        return {
            'wt': len(poss[mask]),
            'del': len(poss_del[mask_del]),
            'overlap': len(overlap_del),
            'overlap_ratio': len(overlap_del) / max(len(poss_del[mask_del]), len(poss[mask])),
            'gain': len(gain),
            'gain_ratio': len(gain) / max(len(poss_del[mask_del]), len(poss[mask])),
            'loss': len(loss),
            'loss_ratio': len(loss) / max(len(poss_del[mask_del]), len(poss[mask]))}


############## Run all scoring functions ###################

def run_scoring_functions(map_a, map_b):
    """
    :param map_a: contact matrix (numpy)
    :param map_b: contact matrix (numpy)
    :return: dictionary of scoring results
    """
    # Basic methods
    mse_score = mse(map_a, map_b)
    spearmanr_score = spearman(map_a, map_b)
    pearsonr_score = pearson(map_a, map_b)
    ssim_score = ssim_map(map_a, map_b)
    scc_score = scc(map_a, map_b)

    # Map-motivated methods
    insulation = vectorMethodToScalar(insulation_track, map_a, map_b)
    contact_decay = vectorMethodToScalar(decay_track, map_a, map_b)
    pca = vectorMethodToScalar(contact_pca_track, map_a, map_b)
    di = vectorMethodToScalar(DI_track, map_a, map_b)

    downres_a = downres(map_a, new_resolution=4681)  # 224 x 224
    downres_b = downres(map_b, new_resolution=4681)  # 224 x 224
    triangle = vectorMethodToScalar(triangle_track, downres_a, downres_b)

    # Biologically motivated methods
    TADs_info = TAD(map_a, map_b, plot=False, detail=True)
    loops_info = loops_diff(map_a, map_b, p=2, width=5, ther=1.1, ther_H=1.1, ther_V=1.1, detail=True)

    return {'mse': mse_score,
            'spearman': spearmanr_score,
            'pearson': pearsonr_score,
            'ssi': ssim_score,
            'scc': scc_score,
            'eigenvector_mse': pca['mse'],
            'eigenvector_spearmanr': pca['spearmanr'],
            'eigenvector_pearsonr': pca['pearsonr'],
            'DI_mse': di['mse'],
            'DI_spearmanr': di['spearmanr'],
            'DI_pearsonr': di['pearsonr'],
            'insulation_mse': insulation['mse'],
            'insulation_spearmanr': insulation['spearmanr'],
            'insulation_pearsonr': insulation['pearsonr'],
            'contact_decay_mse': contact_decay['mse'],
            'contact_decay_spearmanr': contact_decay['spearmanr'],
            'contact_decay_pearsonr': contact_decay['pearsonr'],
            'triangle_mse': triangle['mse'],
            'triangle_spearmanr': triangle['spearmanr'],
            'triangle_pearsonr': triangle['pearsonr'],
            'loops': loops_info,
            'TADs': TADs_info
            }
