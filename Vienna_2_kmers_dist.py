#!/usr/bin/python

## adjusted for data preprocesing (unseen data set for classifier)

"""from (Vienna Format):
>piRNA6
GGGGGAGUCGUUCCAGUAGAUGGAAAG
..........(((((.....))))).. ( -4.50)

#CUGUACAGCCUCCUAGCUUUCC
#......(((......))).... 

to:     
>piRNA6
ctgtacAGCctcctaGCTttcc

final:
csv file (size= total amount of sequences X 8**k_len)+1 :
    get_kmer_dist('ctgtacAGCctcctaGCTttcc')+ lable
"""

import pandas as pd
from for_reut import get_kmer_dist
import numpy as np
import argparse
parser = argparse.ArgumentParser(description="Vienna to RNA bond 8-letter sequence to csv k-mers distribution")
parser.add_argument("-k", default=3, metavar="<value>", type=int, help="kmers length.")
parser.add_argument('-i','--list', action='append', help="Input fold file. - if input is stdin, for train sets positive file first.",required=True) #like: Vienna_2_kmers_dist.py -k 3 -i final_piRNA.fold -i final_nopiRNA.fold
args = parser.parse_args()

k_len = args.k
input_vienna_file = args.list

def eight_letter(rna_seq, dot_bracket):
    """from RNA sequence and dot-bracket notation, get 8-letter representative sequence"""
    rna_out=""
    rna_seq= rna_seq.replace("U", "T")
    for i in range(len(rna_seq)):
        if dot_bracket[i]==".":
            rna_out = rna_out + rna_seq[i].lower()
        else:
            rna_out = rna_out + rna_seq[i]
    return rna_out


def vienna2_arr( file, lable, k_len):
    kmers_B = ''
    print( file , lable)
    with open( file ,"r+") as fold:
        for line in fold:
            if line.startswith(">"):
                #ID = line[1:].rstrip("\n")
                seq = next(fold).rstrip("\n")
                notation = next(fold).rstrip("\n")
                dist = get_kmer_dist(eight_letter(seq, notation), k_len)
                kmers_A = np.array([dist + [lable]]) # bool lables for SL models
                #index_list.append(ID)
                if type(kmers_B)!= str:
                    kmers_B = np.vstack((kmers_B, kmers_A))
                else: 
                    kmers_B = kmers_A.copy() 
    return kmers_B


if len(input_vienna_file) == 1:
    kmers_all = vienna2_arr(input_vienna_file[0], 0, k_len)
    print(np.unique(kmers_all[:,-1], return_counts=True))
    np.savetxt("k%s-mers_dist.csv" %k_len, kmers_all[:,:-1] , delimiter = ",") #without lable. for unseen test set
else:
    kmers_N = vienna2_arr(input_vienna_file[1], 0, k_len) # negative seq
    kmers_P = vienna2_arr(input_vienna_file[0], 1, k_len) #true seq
    kmers_all = np.vstack((kmers_P, kmers_N))
    print(np.unique(kmers_all[:,-1], return_counts=True))
    np.savetxt("k%s-mers_dist_all.csv"%k_len , kmers_all, delimiter = ",")

# =============================================================================
# #                 for feature importance in later models
# from itertools import product
# kmers_count = {''.join(kmer): 0 for kmer in product(['a', 't', 'c', 'g', 'A', 'T', 'C', 'G'], repeat=k_len)}
# feature_names = list(kmers_count.keys())
# =============================================================================

