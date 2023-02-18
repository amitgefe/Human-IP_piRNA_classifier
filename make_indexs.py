# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 15:01:07 2022

@author: Amit
"""

# make an id catalog for merging samples:
import pandas as pd
import os

in_dir = "putative/gene_subtract/"
samples_ls = [s for s in os.listdir(in_dir) if s.endswith(".bed")]

out = ""
for file in samples_ls:
    sample = pd.read_csv(in_dir+file, sep='\t', header=None, 
                 names=['chr', 'start', 'end', 'info', 'XP', 'strand'])
    sample['sampleID'] = file[:-4]
    if len(out)<1:
        out = sample.copy()
    else:
        out = pd.concat([out, sample], ignore_index=True)

    
out = out.reset_index().rename(columns={'index': 'IDS'})
out['IDS'] = out.apply(lambda i : "MR"+str(i.IDS), axis =1)

out.to_csv("putative_ids_index.csv", index=False)
out[['chr', 'start', 'end','IDS','XP', 'strand' ]].to_csv('all4merge.bed', header=False, index=False, sep='\t')

