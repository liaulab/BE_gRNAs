# Created 17 Jan 2023
# Last Edit 14 Apr 2023
# Author: Calvin XiaoYang Hu

import os
import re
import pandas as pd
import numpy as np

DNA_AA_map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
              "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
              "TAT":"Y", "TAC":"Y", "TAA":".", "TAG":".",
              "TGT":"C", "TGC":"C", "TGA":".", "TGG":"W",
              "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
              "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
              "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
              "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
              "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
              "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
              "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
              "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
              "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
              "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
              "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
              "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", }

base_editing_key = {"CBE": ["C", "T"], 
                      "ABE": ["A", "G"]
                     }

bases = 'ACGT'
complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 
                   'a':'t', 't':'a', 'g':'c', 'c':'g'}
cas_key = {'Sp': 'NGG', 'SpG': 'NGN', 'SpRY': 'NNN'}

def DNA_to_AA(seq): 
    aa_seq = ''
    for i in range(len(seq)//3): 
        codon = seq[(i*3):(i*3)+3].replace("T", "U")
        aa_seq += DNA_AA_map[codon]
    return aa_seq

def rev_complement(complements, seq): 
    compl = ''
    for i in range(len(seq)): 
        compl += complements[seq[i]]
    return compl[::-1]

def complement(complements, seq): 
    compl = ''
    for i in range(len(seq)): 
        compl += complements[seq[i]]
    return compl
