# Created 17 Jan 2023
# Author: Calvin XaoYang Hu

import os
import pandas as pd

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
