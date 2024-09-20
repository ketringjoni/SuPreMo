#!/usr/bin/env python
# coding: utf-8


'''
Functions for getting regions of interest (roi).

'''

import pandas as pd
import numpy as np
from pybedtools import BedTool

pred_len = None
d = None

from pathlib import Path
repo_path = Path(__file__).parents[1]

def get_roi(roi, genome):

    if roi == 'genes':
        
        roi_coords = get_TSS_regions(genome, flanking = 1000)
     
    
    else:
        roi_file = pd.read_csv(roi, sep = '\t')
        roi_file['roi_id'] = ['roi_' + str(x) for x in roi_file.index]
        roi_cols = roi_file.columns
        paired_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
        single_cols = ['chr', 'start', 'end']

        if all([x in roi_cols for x in single_cols]):

            roi_coords = roi_file[single_cols + ['roi_id']]

        elif all([x in roi_cols for x in paired_cols]):
    
            paired_loops = roi_file[paired_cols + ['roi_id']]

            # Remove intrachromosomal interactions
            paired_loops_intra = (paired_loops.query('chr1 == chr2')
                                  .rename(columns = {'chr1':'chr'})
                                  .drop('chr2', axis = 1, inplace = False))

            # Remove loops that don't fit in the predictive window 
            paired_loops_intra = remove_large_loops(paired_loops_intra)

            # Group similar loops where the anchors are within distance d
            paired_loops_intra_fit_grouped = group_similar_loops(paired_loops_intra)

            roi_coords = remove_pairing(paired_loops_intra_fit_grouped)

        else:
            raise ValueError('roi file header not compatible.')
    
    return BedTool.from_dataframe(roi_coords.astype({'start': int, 'end': int}))
        

def label_TSS(df):
    if df['strand'] == '+' :
        return df['start']
    if df['strand'] == '-' :
        return df['end']



def get_TSS_regions(genome, flanking = 1000):

    '''
    Get the region around the transcription start site (TSS) of each gene with flanking regions on either end.
    '''

    # Get gene coordinates, protein coding for hg38 and all for hg19 (remove ones with versions, ex AC233755.2)
    gene_annot = pd.read_csv(f'{repo_path}/data/gene_annot_{genome}', sep = ',')
    gene_annot = gene_annot[['.' not in x for x in gene_annot.gene]]

    # Get coordinates for flanking region around TSS
    gene_annot.loc[:,'Start'] = gene_annot.apply(lambda row: label_TSS(row), axis=1) - flanking
    gene_annot.loc[gene_annot['Start'] < 1, 'Start'] = 1
    gene_annot.loc[:,'End'] = gene_annot['Start'] + flanking
    
    return gene_annot[['chr', 'Start', 'End', 'gene']].rename(columns = {'Start':'start', 'End':'end'})





def remove_pairing(loops):

    '''
    Put left and right anchors of paired regions into the same list of regions.
    
    '''
    
    left_anchors = loops[['chr', 'start1', 'end1']].rename(columns = {'start1':'start', 'end1':'end'})
    left_anchors['roi_id'] = loops['roi_id'] + '_L'
    
    right_anchors = loops[['chr', 'start2', 'end2', 'roi_id']].rename(columns = {'start2':'start', 'end2':'end'})
    right_anchors['roi_id'] = loops['roi_id'] + '_R'
    
    return pd.concat([left_anchors, right_anchors], axis = 0).reset_index(drop = True, inplace = False)


    

def remove_large_loops(loops):

    '''
    Remove loops outside of predictive window, using the center of each anchor to measure distance.
    
    '''
    
    center_coord1 = [(x+y)/2 for x,y in zip(loops.start1, loops.end1)]
    center_coord2 = [(x+y)/2 for x,y in zip(loops.start2, loops.end2)]
    
    distance = [abs(x-y) for x,y in zip(center_coord1, center_coord2)]
    
    return loops[[x < pred_len for x in distance]]

    

    

def group_similar_loops(loops):

    '''
    Group loops whose anchors are within d.
    
    '''
    
    left_anchors = loops[['chr', 'start1', 'end1']].drop_duplicates().reset_index(drop = True)
    left_anchors['left_index'] = left_anchors.index
    loops = loops.merge(left_anchors, on = ['chr', 'start1', 'end1'])
    
    right_anchors = loops[['chr', 'start2', 'end2']].drop_duplicates().reset_index(drop = True)
    right_anchors['right_index'] = right_anchors.index
    loops = loops.merge(right_anchors, on = ['chr', 'start2', 'end2'])


    # For all loops with the same left anchor, combine right anchors that are within d
    
    for i in left_anchors.left_index:
        
        loops_i = loops.query('left_index == @i')[['chr', 'start2', 'end2', 'right_index']]
        right_anchors_BED = BedTool.from_dataframe(loops_i)
    
        if len(right_anchors_BED) > 1:
            right_anchor_keep = (right_anchors_BED
                               .sort()
                               .merge(d = d, c = 4, o = 'collapse') # collapse indexes of the anchors merged
                               .to_dataframe()
                               .rename(columns = {'name':'right_index'}))
            
            if len(right_anchor_keep) < len(loops_i):
                
                right_anchor_keep['right_index'] = right_anchor_keep['right_index'].str.split(',')
                right_anchor_keep = right_anchor_keep.explode('right_index')
    
                for ii in right_anchor_keep.right_index.unique():
    
                    loops.loc[(loops.left_index == i) & 
                              (loops.right_index == int(ii)),
                              'start2b'] = int(right_anchor_keep.loc[right_anchor_keep.right_index == ii,'start'])
                    loops.loc[(loops.left_index == i) & 
                              (loops.right_index == int(ii)),
                              'end2b'] = int(right_anchor_keep.loc[right_anchor_keep.right_index == ii,'end'])
                    
                    
    # Repeat for the right anchor
    
    for i in right_anchors.right_index:
        
        loops_i = loops.query('right_index == @i')[['chr', 'start1', 'end1', 'left_index']]
        left_anchors_BED = BedTool.from_dataframe(loops_i)
    
        if len(left_anchors_BED) > 1:
            left_anchor_keep = (left_anchors_BED
                               .sort()
                               .merge(d = 10000, c = 4, o = 'collapse')
                               .to_dataframe()
                               .rename(columns = {'name':'left_index'}))
            
            if len(left_anchor_keep) < len(loops_i):
                
                left_anchor_keep['left_index'] = left_anchor_keep['left_index'].str.split(',')
                left_anchor_keep = left_anchor_keep.explode('left_index')
    
                for ii in left_anchor_keep.left_index.unique():
    
                    loops.loc[(loops.right_index == i) & 
                              (loops.left_index == int(ii)),
                              'start1b'] = int(left_anchor_keep.loc[left_anchor_keep.left_index == ii,'start'])
                    loops.loc[(loops.right_index == i) & 
                              (loops.left_index == int(ii)),
                              'end1b'] = int(left_anchor_keep.loc[left_anchor_keep.left_index == ii,'end'])


    # Fill in the coordinates of anchors that don't change
    for col_name in ['start1b', 'end1b', 'start2b', 'end2b']:
        loops.loc[np.isnan(loops[col_name[:-1]]), col_name] = loops.loc[np.isnan(loops[col_name[:-1]]), col_name[:-1]]


    # Only keep new loop anchor coordinates
    loops_filtered = (loops
                      [['chr', 'start1b', 'end1b', 'start2b', 'end2b', 'roi_id']]
                      .rename(columns = {'start1b':'start1', 'end1b':'end1',
                                         'start2b':'start2', 'end2b':'end2'})
                      .drop_duplicates())


    return remove_large_loops(loops_filtered)









