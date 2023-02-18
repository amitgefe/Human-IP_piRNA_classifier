#!/bin/bash

list=( $(ls | grep fastq.gz | xargs -l basename -s .fastq.gz ) )
## previous: (1185_Cer_L8.LB5 1185_PFX_L7.LB2 1349_Cer_L7.LB3 )   #for pipe test

###                          changed for train set                          ### & pipeline4.sh 
sample_dir="/gpfs0/alal/users/amitgefe/preprocess_chm13v2.0/" #all 266 samples with new genome
##"/gpfs0/alal/users/amitgefe/preprocess/" #all 266 samples
##"/gpfs0/alal/users/amitgefe/preprocess_chm13v2.0/" #all 266 samples with new genome
##"/gpfs0/alal/users/amitgefe/am_pipe/testpip/"   #for pipe test
##"/gpfs0/alal/users/amitgefe/IP_classifier/train_chm13v2.0/" #IP data with new genome
## previous: pipeline1.sh / pipeline3.sh
cd $sample_dir


mkdir -p collapse
mkdir -p hg38_sam
mkdir -p logs/bowtie1
mkdir -p bed/filtered
mkdir -p bed/intersect_ncrna
mkdir -p bed/ncrna_subtract
mkdir -p putative/gene_subtract
mkdir -p logs/split_log
cd logs/split_log

counter=0
for prefix in "${list[@]}"
do
  if [[ "$counter" == '0' ]]; then
    one=$prefix
    ((counter++))
    qsub -cwd -b y -q alal.q -pe ompi 8 -N "st"$one -v prefix=$one /gpfs0/alal/users/amitgefe/am_pipe/supp/pipeline4.sh
    continue
  elif [[ "$counter" == '1' ]]; then
    two=$prefix
    ((counter++))
    qsub -cwd -b y -q alal.q -pe ompi 8 -N "nd"$two -hold_jid "st"$one -v prefix=$two /gpfs0/alal/users/amitgefe/am_pipe/supp/pipeline4.sh
    continue
  fi
    three=$prefix
    qsub -cwd -b y -q alal.q -pe ompi 8 -N "rd"$three -hold_jid "nd"$two -v prefix=$three /gpfs0/alal/users/amitgefe/am_pipe/supp/pipeline4.sh 
  counter=0
done

