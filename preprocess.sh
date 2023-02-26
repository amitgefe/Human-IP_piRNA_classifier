#!/bin/bash
## call in the samples files directory, gz files of fastq.

list=( $(ls | grep fastq.gz | xargs -l basename -s .fastq.gz ) )
mkdir -p preprocess
mkdir -p preprocess/collapse
mkdir -p preprocess/hg38_sam
mkdir -p preprocess/logs/bowtie1
mkdir -p preprocess/bed/filtered
mkdir -p preprocess/bed/intersect_ncrna
mkdir -p preprocess/bed/ncrna_subtract
mkdir -p preprocess/putative/gene_subtract
mkdir -p preprocess/logs/split_log
cd preprocess/logs/split_log

counter=0
for prefix in "${list[@]}"
do
  if [[ "$counter" == '0' ]]; then
    one=$prefix
    ((counter++))
    qsub -cwd -b y -q alal.q -pe ompi 8 -N "st"$one -v prefix=$one /supp/pipeline.sh
    continue
  elif [[ "$counter" == '1' ]]; then
    two=$prefix
    ((counter++))
    qsub -cwd -b y -q alal.q -pe ompi 8 -N "nd"$two -hold_jid "st"$one -v prefix=$two /supp/pipeline.sh
    continue
  fi
    three=$prefix
    qsub -cwd -b y -q alal.q -pe ompi 8 -N "rd"$three -hold_jid "nd"$two -v prefix=$three /supp/pipeline.sh 
  counter=0
done

