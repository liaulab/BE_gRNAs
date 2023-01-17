# Created 17 Jan 2023
# Author: Calvin XaoYang Hu

import os
import pandas as pd

from helper_functions import rev_complement, complement


class BE_gRNAs():

    # references
    bases = 'ACGT'
    complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
                   'a':'t', 't':'a', 'g':'c', 'c':'g'}
    cas_key = {'Sp': 'NGG', 'SpG': 'NGN', 'SpRY': 'NNN'}

    def __init__(self, be_type, editing_window, gene, end_goals, cas_type, exon_filename, exon_dir=''):

        # vars
        self.be_type = be_type
        self.gene = gene
        self.end_goals = end_goals
        self.PAM = self.cas_key[cas_type]
        self.window = editing_window

        # load file of exons
        self.filename = os.path.join(exon_dir, exon_filename)
        f = open(self.filename, "r")
        self.file_content = f.read()
        # list of exons and list of exons with +-20 bps of introns
            # both are necessary bc +-20 bps needed bc gRNA can attach onto DNA that isnt expressed
            # only +- 20 is needed bc gRNA is only 23 bps
        self.exons_extra, self.exons = self.parse_exons()

        if len(self.end_goals) > 0:
            # find all sense and antisense codons that can be mutated
            self.target_codons, self.target_codons_compls = self.generate_target_codons()
            # find all full guide RNAs
            self.gRNAs, self.readdir, self.exon_num, self.targets, self.base_ind = self.get_gRNAs(mode=False)

        else:
            # find all full guide RNAs without a target codon to change
            self.gRNAs, self.readdir, self.exon_num, self.targets, self.base_ind = self.get_gRNAs()

    ###################################################################################################

    # parse exons with and without intron ends from fasta file format
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

    ###################################################################################################

    # generate target codons based on type of base editor, and what codon we want to generate
    def generate_target_codons(self):
        codons = []
        anticodons = []
        for goal in self.end_goals:
            codons += self.base_edits(goal, 'sense')
            anticodons += self.base_edits(goal, 'anti')
        return list(set(codons)), list(set(anticodons))

    def base_edits(self, goal, mode):
        if self.be_type == 'CBE' and mode == 'sense':
            return self.replace(goal, 'C', 'T')
        elif self.be_type == 'CBE' and mode == 'anti':
            return self.replace(goal, 'G', 'A')
        elif self.be_type == 'ABE' and mode == 'sense':
            return self.replace(goal, 'A', 'G')
        elif self.be_type == 'ABE' and mode == 'anti':
            return self.replace(goal, 'T', 'C')
        else: print('BE type invalid')

    # how to do this smarter
    def replace(self, goal, x, y):
        result = [goal[0:i] + x + goal[i+1:] for i in range(3) if goal[i]==y]
        result += [goal[0:i] + x + x + goal[i+2:] for i in range(2) if goal[i:i+2]==y+y]
        result += [x + goal[1] + x for i in range(1) if goal[0]==y and goal[2]==y]
        result += [x+x+x for i in range(1) if goal==y+y+y]
        return result

    ###################################################################################################

    def get_gRNAs(self, mode=True):
        gRNAs = []
        potential_target = []
        readdir = []
        exon_num = []
        base_ind = []
        real_PAM = self.PAM.replace("N", "")

        # make list of all guide RNAs possible
        for i in range(len(self.exons)):
            exon_i = self.exons_extra[i]
            for j in range(len(exon_i)):
                frame = exon_i[j:j+23]
                revcompl_frame = rev_complement(self.complements, frame)

                # looking fwd
                # NGG and NGN cases can be generalized since 21 starts at GG or GN
                if frame[21:21+len(real_PAM)].upper() == real_PAM:
                    sub_frame = frame[self.window[0]-1:self.window[1]]
                    if (self.be_type == 'CBE' and 'C' in sub_frame) or (self.be_type == 'ABE' and 'A' in sub_frame):
                        # no target
                        if mode:
                            gRNAs.append(frame)
                            potential_target.append(frame)
                            readdir.append('fwd')
                            exon_num.append(i+1)
                            x = (len(''.join(self.exons[:i]))) + (j-20) + (frame[3:].find(self.be_type[0])+4) # index of first C
                            base_ind.append(x)
                        # target
                        elif self.checkframe_fwd(frame, i, j) != -1:
                            gRNAs.append(frame)
                            potential_target.append(frame)
                            readdir.append('fwd')
                            exon_num.append(i+1)
                            base_ind.append(self.checkframe_fwd(frame, i, j))

                # looking rev
                # since the opposite just involves looking at the opposite strand, same code but just rev_compl
                if revcompl_frame[21:21+len(real_PAM)].upper() == real_PAM:
                    sub_frame = revcompl_frame[self.window[0]-1:self.window[1]]
                    if (self.be_type == 'CBE' and 'C' in sub_frame) or (self.be_type == 'ABE' and 'A' in sub_frame):
                        # no target
                        if mode:
                            gRNAs.append(revcompl_frame)
                            potential_target.append(frame)
                            readdir.append('rev')
                            exon_num.append(i+1)
                            x = (len(''.join(self.exons[:i]))) + (j-20) + (frame[15:].find(self.be_type[0])+16) # index of last C
                            base_ind.append(x)
                        # target
                        elif self.checkframe_rev(frame, i, j) != -1:
                            gRNAs.append(revcompl_frame)
                            potential_target.append(frame)
                            readdir.append('rev')
                            exon_num.append(i+1)
                            base_ind.append(self.checkframe_rev(frame, i, j))

        return gRNAs, readdir, exon_num, potential_target, base_ind

    ###################################################################################################

    def checkframe_fwd(self, frame, i, j):
        for codon in self.target_codons:
            if codon in frame[self.window[0]-1:self.window[1]+2]:
                ind = frame[self.window[0]-1:self.window[1]+2].find(codon)
                result = len(''.join(self.exons[:i])) + (j-20) + (ind+3)
                if result % 3 == 0:
                    return result
        return -1

    def checkframe_rev(self, frame, i, j):
        for codon in self.target_codons_compls:
            if codon in frame[self.window[0]-1+12:self.window[1]+2+12]:
                ind = frame[self.window[0]-1+12:self.window[1]+2+12].find(codon)
                result = len(''.join(self.exons[:i])) + (j-20) + (ind+15)
                if result % 3 == 0:
                    return result
        return -1

    ###################################################################################################

    def save_data(self, filename):
        head = ['gene', 'guide RNA 5>3', 'target site 5>3', 'direction', 'exon', 'base_ind']
        df = pd.DataFrame(zip([self.gene]*len(self.gRNAs),
                              self.gRNAs,
                              self.targets,
                              self.readdir,
                              self.exon_num,
                              self.base_ind
                             ),
                          columns=head)
        df.to_csv(filename, index=False)
