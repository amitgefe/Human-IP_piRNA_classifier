#!/usr/bin/python

import sys
import pandas as pd
in_file = sys.argv[1]
out_file = sys.argv[2]

nc = "tmp_"+in_file[:-3]+ "txt"
with open( nc , "r") as nc_file:
  ncRNA_lis = [line.rstrip("\n") for line in nc_file]

columns=['chr', 'start_f', 'end_f', 'seq_f', '#reads_x', 'strand', 'CH','start_o', 'end_o', 'read_o', 'XP', 'ST']
intrsct = pd.read_csv(in_file , sep='\t', header=None, names=columns)

intrsct = intrsct[['chr', 'start_f', 'end_f', 'seq_f', '#reads_x', 'strand', 'start_o', 'end_o', 'read_o']]
intrsct["read_o"] = intrsct.read_o.str.replace('::PU', '')
intrsct["start_f"]= intrsct.start_f + 1


def mismatch(seq1, seq2):
    if len(seq1)== len(seq2):
        bool_match_list= [seq1[i]==seq2[i] for i in range(len(seq1))]
        return bool_match_list.count(0)
    else:
        return -1
    
    
intrsct["mm"] = intrsct.apply(lambda x: mismatch(x.seq_f, x.read_o), axis=1)
intrsct["loci_diff"] = intrsct.apply(lambda x: len(set(range(x.start_f, x.end_f)) - (set(range(x.start_o, x.end_o)))), axis=1)
intrsct["nCrNA"] = intrsct.apply(lambda x: ncRNA_lis.count(x.read_o) !=0, axis=1)
intrsct.drop(intrsct[intrsct['seq_f'].isin( intrsct[ intrsct['nCrNA'] == True ].seq_f)].index,inplace=True )

intrsct = intrsct.drop_duplicates(subset=['seq_f', 'read_o'])
intrsct = intrsct[['chr', 'start_f', 'end_f', 'seq_f', '#reads_x', 'strand', 'mm', 'loci_diff', 'read_o']]
print(intrsct.shape)

intrsct.to_csv(out_file, index=False, sep='\t')

