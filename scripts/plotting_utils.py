#!/usr/bin/env python
# coding: utf-8


'''
Functions for plotting predicted contact frequency maps, their disruption tracks, and genes in that region.

'''


# # # # # # # # # # # # # # # # # # 
# # # Import packages and data # # #

import pandas as pd
import math
from cooltools.lib.plotting import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pybedtools import BedTool



# Get Akita measurements (specified in get_scores_utils)
bin_size = 2048
target_length_cropped = 448

map_length = bin_size * target_length_cropped




# # # # # # # # # # # # # # # # # # 
# # # Plotting pre-processing # # #



# Get protein coding gene coordinates (remove ones with versions, ex AC233755.2)
gene_annot_PC = pd.read_csv('../data/gene_annot_PC', sep = ',')
gene_annot_PC = gene_annot_PC[['.' not in x for x in gene_annot_PC.gene]]
gene_annot_PC_BED = BedTool.from_dataframe(gene_annot_PC[['chr', 'start', 'end', 'gene']])





def get_genes_in_map(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN):
    
    
    '''
    Get protein coding genes in map region.
    
    '''

    map_region_BED = BedTool.from_dataframe(pd.DataFrame({'CHR' : [CHR],
                                                          'start' : [map_start_coord],
                                                          'end' : [map_start_coord + map_length]}))

    
    genes_in_map_BED = map_region_BED.intersect(gene_annot_PC_BED, wa = True, wb = True)

    if genes_in_map_BED == '':
        genes_in_map = ''

    else:
        genes_in_map = (genes_in_map_BED
                        .to_dataframe()
                        .rename(columns = {'score': 'Start', 'strand':'End', 'thickStart':'Gene'})
                        [['Start', 'End', 'Gene']])
        
        genes_in_map.Start = [math.ceil(x/bin_size) for x in (genes_in_map.Start - map_start_coord)]
        genes_in_map.loc[genes_in_map.Start < 0,'Start'] = 0
        genes_in_map.End = [math.ceil(x/bin_size) for x in (genes_in_map.End - map_start_coord)]
        genes_in_map.loc[genes_in_map.End > target_length_cropped-1,'End'] = target_length_cropped-1
        genes_in_map['width'] = [x-y+1 for x,y in zip(genes_in_map.End, genes_in_map.Start)]
        
        if SVTYPE == 'DUP' and abs(int(SVLEN)) > bin_size/2:
            
            genes_in_map.loc[(genes_in_map.Start >= rel_pos_map) & 
                              (genes_in_map.End > rel_pos_map),
                              ['Start', 'End']] = genes_in_map.loc[(genes_in_map.Start > rel_pos_map) & 
                                                                  (genes_in_map.End > rel_pos_map),
                                                                  ['Start', 'End']] + math.ceil(SVLEN/bin_size)


            middle_genes = genes_in_map.loc[(genes_in_map.Start < rel_pos_map) & 
                              (genes_in_map.End > rel_pos_map)].index


            for middle_gene in middle_genes:

                left = genes_in_map.loc[[middle_gene]]
                left.End = rel_pos_map
                right = genes_in_map.loc[[middle_gene]]
                right.Start = rel_pos_map + math.ceil(SVLEN/bin_size)
                right.End = right.End + math.ceil(SVLEN/bin_size)

                genes_in_map = (genes_in_map
                                .drop(middle_gene)
                                .append(left)
                                .append(right)
                                .reset_index(drop=True))
            
            # Remove genes that left prediction window
            genes_in_map = genes_in_map[~((genes_in_map.End > target_length_cropped) &
                                        (genes_in_map.Start > target_length_cropped))]
            
            # Crop end of genes that left prediction window
            genes_in_map.loc[(genes_in_map.End > target_length_cropped),'End'] = target_length_cropped
            
            
    return genes_in_map





def get_var_bins(rel_pos_map, SVTYPE, SVLEN):
    
    end = rel_pos_map + math.ceil(SVLEN/bin_size)
    
    if SVTYPE == 'DUP':
        return [rel_pos_map, end, end + math.ceil(SVLEN/bin_size)]
        
    else:    
        return [rel_pos_map, end]
    
    
    
    

# # # # # # # # # # # # # # # # # # 
# # # Plotting functions # # # # #




def plot_maps_genes(REF_pred, ALT_pred, genes_in_map, lines):

    
    gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    
    vmin = -2
    vmax = 2
    linestyle = 'dashed'
    linewidth = 0.3
    color = 'black'

    plt.rcParams['font.size']= 10
    plot_width = 3
    plot_width1D = 1
    fig, gs = gridspec_inches([plot_width], [plot_width, .15, plot_width, .15, plot_width1D])


    plt.subplot(gs[0,0])
    plt.matshow(REF_pred, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.ylabel('Reference matrix',rotation=90)
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color=color, linestyle=linestyle, linewidth = linewidth)
            plt.axhline(y=line, color=color, linestyle=linestyle, linewidth = linewidth)
    plt.gca().yaxis.tick_right()
    plt.xticks([])

    plt.subplot(gs[2,0])
    plt.matshow(ALT_pred, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.ylabel('Alternate matrix',rotation=90)
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color=color, linestyle=linestyle, linewidth = linewidth)
            plt.axhline(y=line, color=color, linestyle=linestyle, linewidth = linewidth)
    plt.gca().yaxis.tick_right()
    plt.xticks([])

    ax1 = fig.add_subplot(gs[4,0])
            
    ax1.broken_barh(gene_track[::3], (14.5, 0.5), facecolors='tab:blue')
    ax1.broken_barh(gene_track[1:][::3], (9.5, 0.5), facecolors='tab:blue')
    ax1.broken_barh(gene_track[2:][::3], (4.5, 0.5), facecolors='tab:blue')
    
    bar_locations = ([14.5, 9.5, 4.5]*math.ceil(len(genes_in_map)/3))[:len(genes_in_map)]
    for i in range(len(genes_in_map)):
        
        bar_location = bar_locations[i]
        
        gene = genes_in_map.loc[i,'Gene']
        location = genes_in_map.loc[i,'Start']
        
        ax1.annotate(gene, (location,bar_location), # annotate gene at the start
                     rotation = 45, 
                     ha = "right", va = "top", # set horizontal and vertical alignment
                     annotation_clip = False, # keep annotation if outside of window
                     fontsize = 8) 
            
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color=color, linestyle=linestyle, linewidth = linewidth)
    plt.ylabel('Genes',rotation=90)
    plt.ylim([0,15])
    plt.xlim([0,448])
    plt.axis('off')

    plt.show()
    
    
    
    
def plot_disruption_tracks(disruption_track, scoring_method):

    plt.figure(figsize=(6,2))

    plt.plot(list(range(target_length_cropped)), disruption_track)
    plt.xlabel('Bins')
    plt.ylabel(scoring_method, rotation = 90)
    plt.xlim([0,target_length_cropped])

    plt.show()
    
    
    
    
def plot_maps_genes_tracks(REF_pred, ALT_pred, genes_in_map, lines, disruption_track, scoring_method):

    gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    
    vmin = -2
    vmax = 2
    linestyle = 'dashed'
    linewidth = 0.3
    color = 'black'

    plt.rcParams['font.size']= 10
    plot_width = 3
    plot_width1D = 1
    fig, gs = gridspec_inches([plot_width], [plot_width, .15, plot_width, .15, 
                                             plot_width1D, .15, plot_width1D*2])


    plt.subplot(gs[0,0])
    plt.matshow(REF_pred, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.ylabel('Reference matrix',rotation=90)
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color=color, linestyle=linestyle, linewidth = linewidth)
            plt.axhline(y=line, color=color, linestyle=linestyle, linewidth = linewidth)
    plt.gca().yaxis.tick_right()
    plt.xticks([])

    plt.subplot(gs[2,0])
    plt.matshow(ALT_pred, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.ylabel('Alternate matrix',rotation=90)
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color=color, linestyle=linestyle, linewidth = linewidth)
            plt.axhline(y=line, color=color, linestyle=linestyle, linewidth = linewidth)
    plt.gca().yaxis.tick_right()
    plt.xticks([])
    
    plt.subplot(gs[4,0])
    plt.plot(list(range(target_length_cropped)), disruption_track)
    plt.ylabel(scoring_method, rotation = 90)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    plt.xticks([])
    

    ax1 = fig.add_subplot(gs[6,0])
            
    ax1.broken_barh(gene_track[::3], (14.5, 0.5), facecolors='tab:blue')
    ax1.broken_barh(gene_track[1:][::3], (9.5, 0.5), facecolors='tab:blue')
    ax1.broken_barh(gene_track[2:][::3], (4.5, 0.5), facecolors='tab:blue')
    
    bar_locations = ([14.5, 9.5, 4.5]*math.ceil(len(genes_in_map)/3))[:len(genes_in_map)]
    for i in range(len(genes_in_map)):
        
        bar_location = bar_locations[i]
        
        gene = genes_in_map.loc[i,'Gene']
        location = genes_in_map.loc[i,'Start']
        
        ax1.annotate(gene, (location,bar_location), # annotate gene at the start
                     rotation = 45, 
                     ha = "right", va = "top", # set horizontal and vertical alignment
                     annotation_clip = False, # keep annotation if outside of window
                     fontsize = 8) 
            
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color=color, linestyle=linestyle, linewidth = linewidth)
    plt.ylabel('Genes',rotation=90)
    plt.ylim([0,15])
    plt.xlim([0,448])
    plt.axis('off')

    plt.show()
    
    
    