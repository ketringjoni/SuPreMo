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
import itertools



# Get Akita measurements (specified in get_Akita_scores_utils)
bin_size = 2048
target_length_cropped = 448

map_length = bin_size * target_length_cropped



# # # # # # # # # # # # # # # # # # 
# # # Plotting pre-processing # # #




# Get gene coordinates, protein coding for hg38 and all for hg19 (remove ones with versions, ex AC233755.2)
gene_annot_hg38 = pd.read_csv('../data/gene_annot_hg38', sep = ',')
gene_annot_hg38 = gene_annot_hg38[['.' not in x for x in gene_annot_hg38.gene]]
gene_annot_hg38_BED = BedTool.from_dataframe(gene_annot_hg38[['chr', 'start', 'end', 'gene']])

gene_annot_hg19 = pd.read_csv('../data/gene_annot_hg19', sep = ',')
gene_annot_hg19 = gene_annot_hg19[['.' not in x for x in gene_annot_hg19.gene]]
gene_annot_hg19_BED = BedTool.from_dataframe(gene_annot_hg19[['chr', 'start', 'end', 'gene']])




def get_genes_in_map(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN, genome = 'hg38'):
    
    
    '''
    Get protein coding genes in map region.
    
    '''

    map_region_BED = BedTool.from_dataframe(pd.DataFrame({'CHR' : [CHR],
                                                          'start' : [map_start_coord],
                                                          'end' : [map_start_coord + map_length]}))

    if genome == 'hg38':
        gene_annot_BED = gene_annot_hg38_BED
    elif genome == 'hg19':
        gene_annot_BED = gene_annot_hg19_BED
        
    
    genes_in_map_BED = map_region_BED.intersect(gene_annot_BED, wa = True, wb = True)

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
            
            genes_in_map.reset_index(drop=True, inplace = True)
            
            
    return genes_in_map





def get_var_bins(rel_pos_map, SVTYPE, SVLEN):
    
    '''
    Get the first and last bin of the variant on the maps.
    
    '''
    
    if SVTYPE == 'BND':
        return [rel_pos_map]
    
    else:
        SVLEN_bins = math.ceil(abs(SVLEN)/bin_size)
        
        if SVLEN_bins == 1:
            return [rel_pos_map]
    
        if SVTYPE == 'DUP':
            return [rel_pos_map, rel_pos_map + SVLEN_bins, rel_pos_map + 2*SVLEN_bins]
        
        else:    
            return [rel_pos_map, rel_pos_map + SVLEN_bins]
    
    
    
def pcolormesh_45deg(plt, mat, lines, linewidth, *args, **kwargs):
    
    '''
    Plot contact frequency maps as triangles (supplement to plot_maps_genes and plot_maps_genes_tracks).
    
    Parts of the code were adapted from https://github.com/pollardlab/contact_map_scoring/blob/main/code/utils.py.
    
    '''
    n = mat.shape[0]
    
    # Create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    
    # Create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    
    # Plot
    im = plt.pcolormesh(A[:,1].reshape(n+1,n+1),A[:,0].reshape(n+1,n+1),np.flipud(mat), *args, **kwargs)
    
    im.set_rasterized(True)
    if lines is not None:
        for line in lines:
            plt.plot([line,line + (448-line)/2], [0,448-line], color = 'gray', linewidth=linewidth, linestyle = 'dashed')
            plt.plot([line,line/2], [0,line], color = 'gray', linewidth=linewidth, linestyle = 'dashed')
    plt.ylim(0,n)
    plt.xticks([])
    plt.yticks([])
    plt.plot([0, n/2], [0, n], 'k-',linewidth=1)
    plt.plot([n/2, n], [n, 0], 'k-',linewidth=1)
    plt.plot([0, n], [0, 0], 'k-',linewidth=1)
    plt.xlim([0,448])

    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_rasterized(True)

    return im




# # # # # # # # # # # # # # # # # # 
# # # Plotting functions # # # # #

 
    
def plot_disruption_tracks(disruption_track, scoring_method):

    '''
        
    Plot disruption score tracks.
    
    '''
    
    plt.figure(figsize=(6,2))

    plt.plot(list(range(target_length_cropped)), disruption_track)
    plt.xlabel('Bins')
    plt.ylabel(scoring_method, rotation = 90)
    plt.xlim([0,target_length_cropped])

    plt.show()
    
    
    

def plot_maps(maps, lines, scale = 1):

        
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. For duplications, there will be 3 lines marking the two regions that are duplicates. 
    
    '''
    
    
    n = maps[0].shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    
    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        
    plt.subplot(gs[4,0])
    pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
        

    plt.show()

    
    
    
def plot_maps_genes(maps, genes_in_map, lines, scale = 1, gene_rows = 3):

        
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. For duplications, there will be 3 lines marking the two regions that are duplicates. 
    
    Plot genes that match the regions in the reference map.
    
    '''
    
    
    gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    
    n = maps[0].shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    
    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width, .15, 
                                               plot_width1D*2*gene_rows/3])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        
    plt.subplot(gs[4,0])
    pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
        

    ax1 = fig.add_subplot(gs[6,0])
    
    gene_bar_height = 15*gene_rows/3
    bar_locations = []
    
    for bar_level in range(gene_rows):

        bar_location = (bar_level+1)*(gene_bar_height/gene_rows) - 0.5
        ax1.broken_barh(gene_track[bar_level:][::gene_rows], 
                        (bar_location, 0.5*scale), 
                        facecolors='tab:blue')
        bar_locations.append(bar_location)

    
    bar_locations = (bar_locations*math.ceil(len(genes_in_map)/gene_rows))[:len(genes_in_map)]

    for i in range(len(genes_in_map)):
        
        bar_location = bar_locations[i]
        
        gene = genes_in_map.loc[i,'Gene']
        location = genes_in_map.loc[i,'Start']
        
        ax1.annotate(gene, (location,bar_location), # annotate gene at the start
                     rotation = 45, 
                     ha = "right", va = "top", # set horizontal and vertical alignment
                     annotation_clip = False, # keep annotation if outside of window
                     fontsize = 12*scale) 
            
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.ylabel('Genes',rotation=90)
    plt.ylim([0,gene_bar_height])
    plt.xlim([0,448])
    plt.axis('off')

    plt.show()

    
    
    
    
def plot_maps_tracks(maps, lines, disruption_track, scale = 1):
    
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. For duplications, there will be 3 lines marking the two regions that are duplicates. 
    
    Plot disruption score tracks.
    
    scale allows you to scale the size of the figure.

    '''

    n = maps[0].shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    
    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width, .15, 
                                               plot_width1D])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        
    plt.subplot(gs[4,0])
    pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
        

    
    plt.subplot(gs[6,0])
    plt.plot(list(range(target_length_cropped)), disruption_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.yticks([])
    

    plt.show()
    
    
    
    

def plot_maps_genes_tracks(maps, genes_in_map, lines, disruption_track, scale = 1, gene_rows = 3):
    
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. For duplications, there will be 3 lines marking the two regions that are duplicates. 
    
    Plot genes that match the regions in the reference map.
    
    Plot disruption score tracks.
    
    scale allows you to scale the size of the figure.
    
    gene_rows allows you to specify the number of rows to plot genes on.
    
    '''

    gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    
    n = maps[0].shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    
    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width, .15, 
                                               plot_width1D, .15, plot_width1D*2*gene_rows/3])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        
    plt.subplot(gs[4,0])
    pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
        

    
    plt.subplot(gs[6,0])
    plt.plot(list(range(target_length_cropped)), disruption_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.yticks([])
    

    ax1 = fig.add_subplot(gs[8,0])
            
    
    gene_bar_height = 15*gene_rows/3
    bar_locations = []
    
    for bar_level in range(gene_rows):

        bar_location = (bar_level+1)*(gene_bar_height/gene_rows) - 0.5
        ax1.broken_barh(gene_track[bar_level:][::gene_rows], 
                        (bar_location, 0.5*scale), 
                        facecolors='tab:blue')
        bar_locations.append(bar_location)

    
    bar_locations = (bar_locations*math.ceil(len(genes_in_map)/gene_rows))[:len(genes_in_map)]

    for i in range(len(genes_in_map)):
        
        bar_location = bar_locations[i]
        
        gene = genes_in_map.loc[i,'Gene']
        location = genes_in_map.loc[i,'Start']
        
        ax1.annotate(gene, (location,bar_location), # annotate gene at the start
                     rotation = 45, 
                     ha = "right", va = "top", # set horizontal and vertical alignment
                     annotation_clip = False, # keep annotation if outside of window
                     fontsize = 12*scale) 
            
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.ylabel('Genes',rotation=90)
    plt.ylim([0,gene_bar_height])
    plt.xlim([0,448])
    plt.axis('off')

    plt.show()
    
    
    
    
    
    
def plot_maps_genes_tracks_nonames(maps, genes_in_map, lines, disruption_track, scale = 1):
    
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. For duplications, there will be 3 lines marking the two regions that are duplicates. 
    
    Plot genes that match the regions in the reference map.
    
    Plot disruption score tracks.
    
    '''

    gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    
    n = maps[0].shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    
    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width, .15, 
                                               plot_width1D, .15, plot_width1D/2])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        
    plt.subplot(gs[4,0])
    pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
        

    
    plt.subplot(gs[6,0])
    plt.plot(list(range(target_length_cropped)), disruption_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.yticks([])
    

    ax1 = fig.add_subplot(gs[8,0])
            
    ax1.broken_barh(gene_track[::3], (1.5, 0.5*scale), facecolors='tab:blue')
    ax1.broken_barh(gene_track[1:][::3], (0.75, 0.5*scale), facecolors='tab:blue')
    ax1.broken_barh(gene_track[2:][::3], (0, 0.5*scale), facecolors='tab:blue')

            
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.ylabel('Genes',rotation=90)
    plt.ylim([0,2])
    plt.xlim([0,448])
    plt.axis('off')

    plt.show()
    
    

    
    
