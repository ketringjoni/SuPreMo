#!/usr/bin/env python
# coding: utf-8



'''
Functions that accompany SuPreMo get_seq for scoring variants for generating mutated sequences.

'''



# # # # # # # # # # # # # # # # # # 
# # # # # Import packages # # # # #

import pandas as pd
import numpy as np

import os
import io

import math
from collections import Counter

from Bio.Seq import Seq



svlen_limit = None
nt = ['A', 'T', 'C', 'G']


seq_length = None
half_patch_size = None




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
    centro_start = int(centromere_coords[centromere_coords.chrom == CHR]['start'])
    centro_stop = int(centromere_coords[centromere_coords.chrom == CHR]['end'])

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
        seq_pos = 0 + adjust + abs(shift) # 1 is added so the position is never 0. coordinates are 1-based
        # POS is the base before variant so 0 is the base before 1
        
    elif position == 'chrom_centro_right':
        seq_pos = centro_stop + adjust + abs(shift)         

    elif position == 'chrom_end':
        seq_pos = chrom_max - seq_length + adjust - abs(shift)       

    elif position == 'chrom_centro_left':
        seq_pos = centro_start - seq_length + adjust - abs(shift)
        
        
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

    REF_half_left = math.ceil((seq_length - REF_len)/2) - shift # if the REF allele is odd, shift right
    REF_half_right = math.floor((seq_length - REF_len)/2) + shift

    
    # Annotate whether variant position with respect to chromosome arms ends
    if len(REF) <= len(ALT):
        var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right)
  
    elif len(REF) > len(ALT):       
        ALT_len = len(ALT)
        ALT_half_left = math.ceil((seq_length - ALT_len)/2) - shift
        ALT_half_right = math.floor((seq_length - ALT_len)/2) + shift   
        var_position = get_variant_position(CHR, POS, ALT_len, ALT_half_left, ALT_half_right)
    

    # Get last coordinate of chromosome
    chrom_max = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max'])
    
    # Get centromere coordinates
    centro_start = int(centromere_coords[centromere_coords.chrom == CHR]['start'])
    centro_stop = int(centromere_coords[centromere_coords.chrom == CHR]['end'])
    
    
    # Get start and end of reference sequence
    if var_position == "chrom_mid":
        REF_start = POS - REF_half_left
        REF_stop = REF_start + seq_length 
    elif var_position == "centromere":
        raise ValueError('Centromeric variant.')
    else:
        REF_start = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, 0, shift)
        REF_stop = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, seq_length, shift)
        print("Warning: Variant not centered; too close to chromosome arm ends.")
        
        
    # Get reference sequence
    REF_seq = fasta_open.fetch(CHR, REF_start, REF_stop).upper()


    # Error if N composition is more than 5% of sequence
    if Counter(REF_seq)['N']/seq_length*100 > 5:
        raise ValueError('N composition greater than 5%.')



    # Variant position relative to the reference sequence outputed
    var_rel_pos_REF = POS - REF_start - 1 
    # subtract 1 to include POS in var_rel_pos_REF since POS is not actually included in the variant

    # Error if reference sequence does not match given REF
    if REF_seq[var_rel_pos_REF : var_rel_pos_REF + REF_len] != REF:
        raise ValueError('Reference allele does not match hg38.')
            
            
            
    # Error if reference sequence is not the right length      
    if len(REF_seq) != seq_length:
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
            ALT_stop = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, seq_length + del_len, shift)
            
        elif var_position in ["chrom_centro_left", "chrom_end"]: 
            ALT_start = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, 0 - del_len, shift)
            ALT_stop = adjust_seq_ends(centro_start, centro_stop, chrom_max, var_position, seq_length, shift)
        
        
        
        # Get alternate sequence
        ALT_seq = fasta_open.fetch(CHR, ALT_start, ALT_stop).upper()
            
            
            
        # Variant position relative to the reference sequence outputed
        var_rel_pos_ALT = POS - ALT_start - 1
        # subtract 1 to include POS in var_rel_pos_REF since POS is not actually included in the variant


        # Error if alternate sequence does not match REF at POS
        if ALT_seq[var_rel_pos_ALT : var_rel_pos_ALT + REF_len] != REF:
            raise ValueError('Sequence for the alternate allele does not match hg38 at REF position.')


    
        # Change alternate sequence to match ALT at POS
        ALT_seq = ALT_seq[:var_rel_pos_ALT] + ALT + ALT_seq[var_rel_pos_ALT + REF_len:]

            
    if len(ALT_seq) != seq_length:
        raise ValueError('Alternate sequence generated is not the right length.')
         
            
    # Take reverse compliment of sequence
    if revcomp:
        REF_seq, ALT_seq = [str(Seq(x).reverse_complement()) for x in [REF_seq, ALT_seq]]

        
    return REF_seq, ALT_seq, [var_rel_pos_REF, var_rel_pos_ALT]










# # # # # # # # # # # # # # # # # # 
# # # Generating BND sequences # # #



def adjust_seq_ends_BND(CHR, position, adjust, shift):
           
    '''
    Get start (adjust = 0) or end (adjust = seq_length) of sequence for prediction based on variant position \
    with respect to chromosome arm ends (defined in get_variant_position function).
    Different from adjust_seq_ends because it does not require centro_start, centro_stop, and chrom_max as input.
    1 is added to convert position to 1-based.
    
    '''
    
    if position == 'chrom_start':
        seq_pos = adjust + abs(shift) + 1 
        
    elif position == 'chrom_centro_right':
        seq_pos = int(centromere_coords[centromere_coords.chrom == CHR]['end']) + adjust + abs(shift) + 1

    elif position == 'chrom_end':
        seq_pos = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max']) - seq_length + adjust - abs(shift) + 1
        
    elif position == 'chrom_centro_left':
        seq_pos = int(centromere_coords[centromere_coords.chrom == CHR]['start']) - seq_length + adjust - abs(shift) + 1
        
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

        difference = adjust_seq_ends_BND(self.CHR_right, self.right_position, seq_length, self.shift) - self.right_end

        self.left_start += difference
        self.right_end += difference

        return self.left_start, self.right_end


    def adjust_both_start(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, seq_length, self.shift)

        difference = (self.POS_left - left_start) - (self.POS_right - (right_end - seq_length))

        if difference <=0:
            # left side is closest to start (chromosome start or right end of centromere)
            right_end -= difference

        else:
            left_start += difference

        return left_start, right_end


    def adjust_both_end(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, seq_length, self.shift)

        difference = (left_start + seq_length - self.POS_left) - (right_end - self.POS_right)

        if difference <=0:
            # left side is closest to end (chromosome end or left end of centromere)
            right_end += difference

        else:
            left_start -= difference

        return left_start, right_end




    # Antisense left


    def adjust_left_antisense_left(self):

        difference = adjust_seq_ends_BND(self.CHR_left, self.left_position, seq_length, self.shift) - self.left_start

        self.left_start += difference 
        self.right_end -= difference

        return self.left_start, self.right_end



    def adjust_right_antisense_left(self):

        difference = adjust_seq_ends_BND(self.CHR_right, self.right_position, seq_length, self.shift) - self.right_end

        self.left_start -= difference
        self.right_end += difference 

        return self.left_start, self.right_end



    def adjust_start_end_antisense_left(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, seq_length, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, seq_length, self.shift)

        difference = (self.POS_left - (left_start - seq_length)) - (right_end - self.POS_right)

        if difference <=0:
            # left side is closest to edge
            right_end += difference

        else:
            left_start += difference

        return left_start, right_end



    def adjust_end_start_antisense_left(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, seq_length, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, seq_length, self.shift)

        difference = (left_start - self.POS_left) - (self.POS_right - (right_end - seq_length))

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

        difference = (self.POS_left - left_start) - (seq_length - (self.POS_right - right_end))

        if difference <=0:
            # left side is closest to edge
            right_end += difference

        else:
            left_start += difference

        return left_start, right_end




    def adjust_end_start_antisense_right(self):

        left_start = adjust_seq_ends_BND(self.CHR_left, self.left_position, 0, self.shift)
        right_end = adjust_seq_ends_BND(self.CHR_right, self.right_position, 0, self.shift)

        difference = (left_start + seq_length - self.POS_left) - (self.POS_right - right_end)

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
    
    if BND_rel_pos != seq_length/2:
        print("Warning: Variant not centered; too close to chromosome arm ends.")
    
    
    # Get sequences with adjusted sequence ends
    ALT_left = fasta_open.fetch(CHR_left, left_start, POS_left).upper()

    ALT_right = fasta_open.fetch(CHR_right, POS_right - 1, right_end - 1).upper() 

    REF_for_left = fasta_open.fetch(CHR_left, left_start, left_start + seq_length).upper()
    REF_for_right = fasta_open.fetch(CHR_right, right_end - 1 - seq_length, right_end - 1).upper() 
          
    
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
    
    if BND_rel_pos != seq_length/2:
        print("Warning: Variant not centered; too close to chromosome arm ends.")
      
        
    # Get sequences with adjusted sequence ends
    ALT_left_revcomp = fasta_open.fetch(CHR_left, POS_left - 1, left_start - 1).upper()
    ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())

    ALT_right = fasta_open.fetch(CHR_right, POS_right - 1, right_end - 1).upper()

    REF_for_left_revcomp = fasta_open.fetch(CHR_left, left_start - 1 - seq_length, left_start - 1).upper() 
    REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
    REF_for_right = fasta_open.fetch(CHR_right, right_end - 1 - seq_length, right_end - 1).upper() 
          
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
    
    if BND_rel_pos != seq_length/2:
        print("Warning: Variant not centered; too close to chromosome arm ends.")
      

    # Get sequences with adjusted sequence ends
    ALT_left = fasta_open.fetch(CHR_left, left_start, POS_left).upper()

    ALT_right_revcomp = fasta_open.fetch(CHR_right, right_end, POS_right).upper()
    ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

    REF_for_left = fasta_open.fetch(CHR_left, left_start, left_start + seq_length).upper()
    REF_for_right_revcomp = fasta_open.fetch(CHR_right, right_end, right_end + seq_length).upper()
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
        raise ValueError('Centromeric variant.')

    if '[' in ALT:
        
        CHR2 = ALT.split(':')[0].split('[')[1]
        POS2 = int(ALT.split('[')[1].split(':')[1])
        
        var_position2 = get_variant_position(CHR2, POS2, 0, half_patch_size - shift, half_patch_size + shift)
        
        if var_position2 == 'centromere':
            raise ValueError('Centromeric variant.')

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
            raise ValueError('Centromeric variant.')

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
            raise ValueError('Unexpected format: BND ALT does not include REF.')
    
        # Create ALT sequence
        ALT_seq = ALT_left + ALT_t + ALT_right
        
        # chop off the sides if longer than seq_length
        ALT_seq = crop_sequence(ALT_seq, seq_length)
        
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
            raise ValueError(f'Variant larger than set limit.')
        
        sequences = get_sequences(CHR, POS, REF, ALT, shift, revcomp)
               
    elif SVTYPE == "BND":

        sequences = get_sequences_BND(CHR, POS, REF, ALT, shift, revcomp)
        
    else:
        raise ValueError('SV type not supported.')
      
    return sequences
















