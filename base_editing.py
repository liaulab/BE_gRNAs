# Created 17 Jan 2023
# Last Edit 14 Apr 2023
# Author: Calvin XiaoYang Hu

import os
import re
import pandas as pd
import numpy as np

from helper_functions import DNA_AA_map, base_editing_key, bases, complements, cas_key
from helper_functions import DNA_to_AA, rev_complement, complement

class Gene_BaseEditing(): 
    
    def __init__(self, filepath, output_dir=''): 
        f = open(filepath, "r")
        self.file_content = f.read()
        self.exons_extra, self.exons = self.parse_exons()
    
    # function to read in a .fasta file as its exons, separating into just exon and exon +- 20 bps
    # outputs: list of exons, list of exons +- 20 bps
    def parse_exons(self): 
        exons_extra = []
        i = -1
        for line in self.file_content.split('\n'): 
            if len(line) > 0 and line[0] == '>': 
                exons_extra.append('')
                i += 1
            else: 
                exons_extra[i] += line
        exons = []
        for exon in exons_extra: 
            exons.append(''.join([base for base in exon if base.isupper()]))
        return exons_extra, exons
    
    # generates all fwd and rev potential guides and their metadata
    # while this is computationally naive, genes are usually not that long
    # pseudo outputs: guide sequence, the index of the first bp, the frame of the first bp, the exon # of guide
    def find_all_guides(self, n): 
        self.fwd_guides = []
        prev_frame = 0
        prev_ind = 0
        for e, exon_extra in enumerate(self.exons_extra): 
            for i in range(len(exon_extra)-22): 
                frame = (i+prev_frame-20)%3
                ind = i-20+prev_ind
                self.fwd_guides.append([exon_extra[i:i+23], frame, ind, e])
            prev_frame = (prev_frame+len(exon_extra)-40)%3
            prev_ind += len(exon_extra)-40
            
        self.rev_guides = [[rev_complement(complements, 
                                           g[0])] + [(g[1]+1)%3] + [g[2]+22] + [g[3]] for g in self.fwd_guides]

        
def find_gRNAs(gene_object, mode, cas_type, target_codons=[], window=[4,8], PAM=None): 
    # mode: can be CBE or ABE, or just a list of 2 base edits ex ["C", "T"]
    # cas_type: Sp, SpG, SpRY
    # target_codons: list of codons that we want to make with our base edit
    # window: 4th to 8th bases inclusive by default, can be changed
    # PAM: optional field to input a custom PAM
    # Returns: a df of exon #, guides (23 bps), target (20 bps), fwd or rev, codon #, edit made
    
    # process mode
    if len(mode) == 3: 
        mode = base_editing_key[mode]
    assert len(mode) == 2
    
    # process PAM
    if PAM is None: 
        PAM = cas_key[cas_type]
    PAM_regex = process_PAM(PAM)
    
    # filter for PAM and contains editable base in window
    fwd_results = [g.copy() for g in gene_object.fwd_guides if PAM_regex.match(g[0][-len(PAM):]) and 
                                                    mode[0] in g[0][window[0]-1:window[1]]]
    for g in fwd_results: 
        # mutates all residues according to the mode
        original = g[0][1:12]
        mutated = g[0][1:window[0]-1] + g[0][window[0]-1:window[1]].replace(mode[0], 
                                                                            mode[1]) + g[0][window[1]:window[1]+4]
        assert(len(original)==len(mutated))
        # compares the residues to find which amino acids were altered and catalogs them
        start = (-1*(g[1]-1))+1
        original, mutated = original[start:start+9], mutated[start:start+9]
        original_aa, mutated_aa = '', ''
        edit, edit_ind = [], []
        for i in range(3): 
            if not original[i*3:(i+1)*3].isupper(): 
                continue
            original_aa += DNA_AA_map[original[i*3:(i+1)*3]]
            mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
            if original_aa[-1] != mutated_aa[-1]: 
                edit.append(original_aa[-1] + ">" + mutated_aa[-1])
                edit_ind.append(int((g[2]+1+start+(i*3))/3))
                
        if len(edit) == 0: 
            edit.append('No Change')
        # append all information to dataframe
        g.append(original_aa)
        g.append(mutated_aa)
        g.append(edit)
        g.append(edit_ind)
        g.append('fwd')
        assert(len(g)) == 9
        
    # filter for PAM and contains editable base in window
    rev_results = [g.copy() for g in gene_object.rev_guides if PAM_regex.match(g[0][-len(PAM):]) and 
                                                    mode[0] in g[0][window[0]-1:window[1]]]
    for g in rev_results: 
        # mutates all residues according to the mode
        original = g[0][1:12]
        mutated = g[0][1:window[0]-1] + g[0][window[0]-1:window[1]].replace(mode[0], 
                                                                            mode[1]) + g[0][window[1]:window[1]+4]
        assert(len(original)==len(mutated))
        # compares the residues to find which amino acids were altered and catalogs them
        original = rev_complement(complements, original[g[1]:g[1]+9])
        mutated = rev_complement(complements, mutated[g[1]:g[1]+9])
        original_aa, mutated_aa = '', ''
        edit, edit_ind = [], []
        for i in range(3): 
            if not original[i*3:(i+1)*3].isupper(): 
                continue
            original_aa += DNA_AA_map[original[i*3:(i+1)*3]]
            mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
            if original_aa[-1] != mutated_aa[-1]: 
                edit.append(original_aa[-1] + ">" + mutated_aa[-1])
                edit_ind.append(int((g[2]-5+(i*3))/3))
                
        if len(edit) == 0: 
            edit.append('No Change')
        # append all information to dataframe
        g.append(original_aa)
        g.append(mutated_aa)
        g.append(edit)
        g.append(edit_ind)
        g.append('rev')
        assert(len(g)) == 9
        
    results = fwd_results + rev_results
    return pd.DataFrame(results)

# function to change a PAM sequence into a regex sequence
def process_PAM(PAM): 
    PAM = PAM.replace("G", "[gG]{1}")
    PAM = PAM.replace("C", "[cC]{1}")
    PAM = PAM.replace("T", "[tT]{1}")
    PAM = PAM.replace("A", "[aA]{1}")
    PAM = PAM.replace("N", "[acgtACGT]{1}")
    return re.compile("({})".format(PAM))
    