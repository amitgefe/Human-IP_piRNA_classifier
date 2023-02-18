##done setpbysetp in /gpfs0/alal/users/amitgefe/preprocess_chm13v2.0/
#!/bin/bash
##      preprocessing test set for XGB classifier      ##

### call from /gpfs0/alal/users/amitgefe/preprocess_chm13v2.0/logs/ with qsub -cwd -b y -q alal.q 
tools_dir="/gpfs0/alal/users/amitgefe/am_pipe/supp/" 
cd ..
out_dir=`pwd` 
echo $out_dir
date

#23/7 redo from putative/gene_subtract/.bed to XGB_result
##merge by genomic loci between samples
/gpfs0/alal/projects/utils/python3.8 $out_dir/make_indexs.py #make one bed file , index file (all4merge.bed, putative_ids_index.csv)
sort -k1,1 -k2,2n all4merge.bed > in_sorted.bed
/gpfs0/alal/projects/utils/bedtools2/bin/bedtools merge -i in_sorted.bed -s -c 4,5,6 -o distinct,sum,distinct > gene_subtract_merged.bed

##get sequences & filter length  . 13,918 seqs. 331,010 loci.
awk 'OFS="\t" { if ($3-$2+1 > 36) {next} ; $2=$2-1} {print $0}' gene_subtract_merged.bed | /gpfs0/alal/projects/utils/bedtools2/bin/bedtools getfasta  -s -name+ -fi "/gpfs0/alal/users/amitgefe/genome_chm13v2/chm13v2.0.genome.fa" -bed -  -fo gene_subtract_merged.fasta

#collapse all seqs
cat gene_subtract_merged.fasta | grep -v ">" | sort -u |  awk '{if (/'--'/) {next}}{print ">seq"NR"|PU\n"$0}' > $out_dir/collapse_putative_updated.fasta
#get fold
/gpfs0/alal/projects/utils/ViennaRNA-2.4.14/bin/RNAfold -i $out_dir/collapse_putative_updated.fasta  > $out_dir/collapse_putative_updated.fold
find ./ -name "*.ps" -delete

#k-mers probability distrubution
/gpfs0/alal/projects/utils/python3.8 $tools_dir/Vienna_2_kmers_dist.py -k 3 -i $out_dir/collapse_putative_updated.fold


/gpfs0/alal/projects/utils/python3.8 $out_dir/make_results_updated.py "/gpfs0/alal/users/amitgefe/IP_classifier/XGB_model.sav" $out_dir
#17M Jul 24 10:58 XGB_results_updated.csv

echo "done" 
date