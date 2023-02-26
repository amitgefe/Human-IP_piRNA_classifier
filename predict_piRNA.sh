#!/bin/bash
##      preprocessing new set for the human piRNA XGB classifier      ##

### call from /preprocess/logs/  
tools_dir="/supp/" 
cd ..
out_dir=`pwd` 
echo $out_dir
date


##merge by genomic loci between samples
python3.8 $tools_dir/make_indexs.py #make one bed file , index file (all4merge.bed, putative_ids_index.csv)
sort -k1,1 -k2,2n all4merge.bed > in_sorted.bed
bedtools merge -i in_sorted.bed -s -c 4,5,6 -o distinct,sum,distinct > gene_subtract_merged.bed

##get sequences & filter length
awk 'OFS="\t" { if ($3-$2+1 > 36) {next} ; $2=$2-1} {print $0}' gene_subtract_merged.bed | bedtools getfasta  -s -name+ -fi "/genome_chm13v2/chm13v2.0.genome.fa" -bed -  -fo gene_subtract_merged.fasta

#collapse all seqs
cat gene_subtract_merged.fasta | grep -v ">" | sort -u |  awk '{if (/'--'/) {next}}{print ">seq"NR"|PU\n"$0}' > $out_dir/collapse_putative_updated.fasta

#get fold
ViennaRNA-2.4.14/bin/RNAfold -i $out_dir/collapse_putative_updated.fasta  > $out_dir/collapse_putative_updated.fold
find ./ -name "*.ps" -delete

#k-mers probability distrubution
python3.8 $tools_dir/Vienna_2_kmers_dist.py -k 3 -i $out_dir/collapse_putative_updated.fold


python3.8 $tools_dir/make_results_updated.py "$tools_dir/XGB_model.sav" $out_dir


echo "done" 
date