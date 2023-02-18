#!/usr/bin/python
## original ~/IP_data/preprocess/putative/preClassifier/generate_conf_piRNAs.py
## copy of "/gpfs0/alal/users/amitgefe/IP_classifier/train/generate_conf_piRNAs.py"

import pandas as pd
import os
import sys

path= os.getcwd()
IP_file = sys.argv[1]
in_dir = sys.argv[2]
#in_dir="/gpfs0/alal/users/amitgefe/IP_classifier/train/putative/gene_subtract/"
# previous - "/gpfs0/alal/users/amitgefe/IP_data/preprocess/putative/gene_subtract/"
os.chdir(in_dir)

IP = pd.read_csv(IP_file, sep='\t', header=None, 
                 names=['chr', 'start', 'end', 'count', '#reads', 'strand'])
                 

ctrl = pd.read_csv("SRR2157868.bed", sep='\t', header=None, 
                   names=['chr', 'start', 'end', 'count', '#reads', 'strand'])
                                                          
IP['source'] = 'IP'
ctrl['source'] = 'ctrl'

IP = IP[(IP.end - IP.start + 1) <= 36] ## changed from 33 to 38. 25/4/22 to 36

ctrl = ctrl[(ctrl.end - ctrl.start + 1) <= 36] ## changed from 33 to 38. 25/4/22 to 36

mrg = pd.merge(IP, ctrl, on=['chr', 'start', 'end', 'strand'], how='outer')                                                               

mrg = mrg.sort_values(by=['chr', 'start', 'end', 'strand'])


intrsct = pd.DataFrame(columns=['chr', 'start', 'end', 'count_x', '#reads_x', 'strand', 'source_x','count_y', '#reads_y', 'source_y'])

i = 0
while i <= mrg.shape[0]:
    first_series = mrg.iloc[i]
    try:
        second_series = mrg.iloc[i+1]
    except IndexError:
        pass
    if (first_series.chr == second_series.chr and 
        first_series.strand == second_series.strand and 
        set(range(first_series['start'], first_series['end']+1)).intersection(set(range(second_series['start'], second_series['end']+1)))):
        chromosome = first_series['chr']
        strand = first_series['strand']
        if first_series['source_x'] == 'IP':
            start = first_series['start']
            end = first_series['end']
            count_x = first_series['count_x']
            reads_x = first_series['#reads_x']
            source_x = first_series['source_x']
            count_y = second_series['count_y']
            reads_y = second_series['#reads_y']
            source_y = second_series['source_y']
        else:
            start = second_series['start']
            end = second_series['end']
            count_x = second_series['count_x']
            reads_x = second_series['#reads_x']
            source_x = second_series['source_x']
            count_y = first_series['count_y']
            reads_y = first_series['#reads_y']
            source_y = first_series['source_y']       
        intrsct = intrsct.append(dict(zip(['chr', 'start', 'end', 'count_x', '#reads_x', 'strand', 'source_x','count_y', '#reads_y', 'source_y'],
                                           [chromosome, start, end, count_x, reads_x, strand, source_x, count_y, reads_y, source_y])), ignore_index=True)
        i += 2
    else:
        intrsct = intrsct.append(mrg.iloc[i])
        i += 1
            

intrsct = intrsct[intrsct['source_x'] == 'IP']

#change all count_x/y to #reads_x/y
intrsct['#reads_x'] = intrsct['#reads_x'].fillna(1e-6)
intrsct['#reads_y'] = intrsct['#reads_y'].fillna(1e-6)

intrsct['norm_count_x'] = (intrsct['#reads_x']/intrsct['#reads_x'].sum())*1e6
intrsct['norm_count_y'] = (intrsct['#reads_y']/intrsct['#reads_y'].sum())*1e6


intrsct['IP_ctrl_ratio'] = intrsct['norm_count_x']/intrsct['norm_count_y']

sig_PI = intrsct[intrsct['IP_ctrl_ratio'] > 1]

sig_PI = sig_PI[['chr', 'start', 'end', 'count_x', 'IP_ctrl_ratio', 'strand']]
sig_PI.columns = ['#chr', 'start', 'end', 'raw_count', 'IP_ctrl_ratio', 'strand']

os.chdir(path)
sig_PI.to_csv(IP_file[:IP_file.find('.bed')] + '.piRNA.bed', index=False, sep='\t')
