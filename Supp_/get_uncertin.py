#!/usr/bin/python
# get all nopiRNA reads minus reads in seqmap ,& piRNA seq in seqmap

import pandas as pd
seqmap = pd.read_csv("seqmap_seq.out" , sep='\t')
with open( "nopiRNA_all.fasta" , "r+") as nopi_file: 
  nopi_dict = {line.rstrip("\n").lstrip(">"): next(nopi_file).rstrip("\n") for line in nopi_file if line.startswith('>')}

seqmap = seqmap[(seqmap.strand == "+") & (seqmap.trans_coord < 3) & (seqmap.num_mismatch < 3)]
seqmap["nopi_seq"] = seqmap.apply(lambda x: nopi_dict[x.trans_id], axis=1)

#pi:
seqmap2 = seqmap[['probe_id','probe_seq' ]].drop_duplicates()
seqmap2["probe_id"]= ">" + seqmap2.probe_id
seqmap2.to_csv("tmp_uncertin_pi.fasta", index=False, header= False, sep='\n')

#nopi:
seqmap1 = seqmap[['trans_id','nopi_seq' ]].drop_duplicates()
seqmap1["trans_id"]= ">" + seqmap1.trans_id
seqmap1.to_csv("tmp_uncertin_nopi.fasta", index=False, header= False, sep='\n')
