#!/usr/bin/python
## make an informative results df from XGB output.

import os
import sys
import numpy as np
import pandas as pd
import pickle
import re


def FASTA_reader(a_file):
    """Reads fasta files to dictionary"""
    dic={}
    with open (a_file,'r') as file:
        for line in file:
            if line.startswith('>'):
                name=line.rstrip('\n').lstrip('>') #read ID
                dic[name]=''            #read ID as key
            else:
                dic[name]+=line.rstrip('\n')      #every line after is part of the sequence
    return dic
   
def find_loci_updated(seq):
    ''' find which samples file contain the output seq 
    get loci from informative fasta header. update 24/7/22
    output per seq (str):
    'chr13\t934483\t344859\t+, chr15\t470424\t704269\t+, chr22\t571738\t717410\t+'
    '''
    loci = os.popen("grep -w -B 1 %s gene_subtract_merged.fasta | grep '>.*)$' | cut -f 3,4 -d ':' -s| sort -u |  tr ':(' '\t'"%seq).read()
    ls = loci.replace("-)", "N").replace("-", "\t").replace("N", "-)").split(")\n")[:-1]
    if len(ls)==0:
        raise ValueError("couldn't find seq in putative/gene_subtract/.fasta")
    return ", ".join(ls) 
    

merged_df = pd.read_csv("putative_ids_index.csv", usecols=['IDS', 'sampleID', 'info'] ) #putative_ids_index.csv
print("merged_df shape " , merged_df.shape)
def find_samples_updated(seq):
    ''' find which samples file contain the output seq . 
        output per seq (str):
        'AN08677_ba41_42_22_1st::XP:3::NR:1, AN17425_ba9_31::XP:3::NR:1'
        update 14/7/22
    '''
    global merged_df 
    header = os.popen("grep -w -B 1 %s gene_subtract_merged.fasta | grep '>.*)$' | cut -f 1 -d ':' -s | tr -s ',>' '\n' | sort -u"%seq).read()
    header = header.split("\n")[1:]
    mask = merged_df.IDS.isin(header)
    dfs = merged_df[mask]
    dfs= dfs[[ 'sampleID','info']].drop_duplicates().apply( lambda s : "::".join(s.values) , axis=1)
    merged_df.drop(mask[mask==True].index , inplace=True)
    samples_ls = dfs.tolist()
    if len(samples_ls)==0:
        raise ValueError("couldn't find seq in putative/gene_subtract/.fasta")
    return ", ".join(samples_ls)


print("\n\n" , os.getcwd())
model_path = sys.argv[1]
# load the model from disk:
loaded_model = pickle.load(open(model_path , 'rb'))
print( "loaded_model.get_params() = \n", loaded_model.get_params())
# predict sequences from kmers csv file:
csv_path = sys.argv[2] #./
in_file = csv_path + "k3-mers_dist.csv"
data = np.genfromtxt( in_file ,delimiter=',' )
print("data shape" , data.shape)
y_pred = loaded_model.predict(data) 

# dictionary of sequences, info header as keys:
id2seq = FASTA_reader(csv_path +  "collapse_putative_updated.fasta") 

# =============================================================================
# # make result df
# =============================================================================
df = pd.DataFrame(y_pred.T , columns=(['pred']))
df['seqID'] = df.reset_index().apply(lambda x: "seq"+ str(x['index']+1) +"|PU", axis =1)
df['seq'] = df.seqID.apply(lambda x: id2seq[x] if x in id2seq else None) 

df['samples'] = df.seq.apply(lambda s: find_samples_updated(s) )
df['loci'] = df.seq.apply(lambda s: find_loci_updated(s) )

df.to_csv( "XGB_results_updated.csv")
print("merged_df shape " , merged_df.shape)
print( "y_pred " , np.unique(y_pred , return_counts=True))
