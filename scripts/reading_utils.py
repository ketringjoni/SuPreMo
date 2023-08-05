#!/usr/bin/env python
# coding: utf-8


'''
Functions for reading input.

'''


import pandas as pd
import numpy as np

import io

import gzip


var_set_size = None





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
                variants['CIEND'] = variants.INFO.str.split('CIEND=').str[1].str.split(';').str[0]
                variants.loc[variants.END == variants.CIEND,'END'] = np.nan
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
        
        colnames = ['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN']
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

            
    elif 'txt' in in_file:
        
        colnames = pd.read_csv(in_file, sep = '\t', nrows = 0).columns
        cols_to_use = ['CHROM', 'POS', 'REF', 'ALT']

        if 'SVTYPE' in colnames:
            cols_to_use.append('END')
            cols_to_use.append('SVTYPE')
            cols_to_use.append('SVLEN')

        variants = pd.read_csv(in_file, sep = '\t', low_memory=False, names = colnames,
                               skiprows = 1 + var_set*var_set_size, nrows = var_set_size)[cols_to_use]
    
    else:
        raise ValueError('Input file type not accepted. Make sure it has the right extension.')
        
        
    variants.reset_index(inplace = True, drop = True)
    
    
    return variants



