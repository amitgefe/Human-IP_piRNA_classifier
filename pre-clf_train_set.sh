#!/bin/bash
##after preprocessing of IP data, make final positive & negative sets.
## original ~/IP_classifier/train/final_VS_origin.sh
#copy of  ~/IP_classifier/train_chm13v2.0/final_VS_origin.sh & ~/IP_classifier/train_chm13v2.0/final_set/final_step.sh

#call from ~/IP_classifier/train_chm13v2.0/logs/ 
cd ..
dir2=`pwd` #"/gpfs0/alal/users/amitgefe/preprocess"
echo $dir2
date
tools_dir="/gpfs0/alal/users/amitgefe/am_pipe/supp/" 
#previous- "/gpfs0/alal/users/amitgefe/IP_classifier/train"
original_dir=$dir2"/collapse/"


###                          changed 24/4 - added | tr -d "::PU"
cat bed/intersect_ncrna/*.txt | tr -d "::PU" > bed/intersect_ncrna/all.ID.txt
gunzip $original_dir*
mkdir final_set
cd final_set
list=(SRR2157869 SRR2157596 SRR2157870)


for prefix in "${list[@]}"
do
  echo $prefix
###                          changed 25/4 - generate_conf_piRNAs.py size selection 36.
  # 1: get final bed with XP::NR:: info :
  /gpfs0/alal/projects/utils/python3.8 $tools_dir/generate_conf_piRNAs.py $prefix.bed $dir2/putative/gene_subtract/

  # 2: get final fa with the informative header
  awk 'OFS="\t" {$2=$2-1; print $0 }' $prefix.piRNA.bed  | /gpfs0/alal/projects/utils/bedtools2/bin/bedtools getfasta  -s -name+  -fi /gpfs0/alal/users/amitgefe/genome_chm13v2/chm13v2.0.genome.fa -bed -  -fo $prefix.piRNA.fasta
  ### 2.5 : list of seq without full match to original reads:
  grep -v ">" $prefix.piRNA.fasta | grep -vwxFf $original_dir$prefix* > $prefix.unmtch2ori.txt
  
  # 3: get seq fasta with the informative header (unmatched to original reads 100%):
  grep -B 1 -wxFf $prefix.unmtch2ori.txt $prefix.piRNA.fasta | awk '{if (/'--'/) {next}}{print $0}'  > $prefix.unmtch2ori.fasta
###                          changed 24/4 - added 3.5.
  ### 3.5. check also seq that are FM and in MR loci.
  grep -B 1 -wxFf $original_dir$prefix* $prefix.piRNA.fasta | tr -s ":()" "\t" | awk ' {if ($4 ==1) {next} {print $0}}' | grep ">" -A 1  | awk '{if (/'--'/) {next}}{print $0}' >> $prefix.unmtch2ori.fasta

###                          changed 25/4 - add +1 to $6.
  # 4: transform all unmatched to bed file & fit to bedtools:
  cat $prefix.unmtch2ori.fasta | tr -s ":()-" "\t" | awk 'OFS = "\t" {{if (NF==7) ST = "-"; else ST = $8}; getline seq ; {print $5, $6+1, $7 ,seq, ($3$4), ST}}' | tr -d ">" > $prefix.fitd.bed
  
  # 5: intersct seq to reads :
  /gpfs0/alal/projects/utils/bedtools2/bin/bedtools intersect -s -wa -wb -nonamecheck -a $prefix.fitd.bed -b  $dir2"/bed/filtered/"$prefix.bed > $prefix.intrsct.bed  
  # 6: tag reads intersect to ncRNAs & explore :
  grep -wFf $dir2/bed/intersect_ncrna/all.ID.txt $prefix.intrsct.bed | cut -f 10 | sort -u | tr -d "::PU" > tmp_$prefix.intrsct.txt 
  /gpfs0/alal/projects/utils/python3.8 $tools_dir/intrsct_explor.py $prefix.intrsct.bed $prefix.exp.bed

  # 7: make reads list ot piRNA (fitted full match & not intersect ncRNAs OR in loci of origin merged seq OR have 1-2 mm & 0 loci diff to seq):
  #7.1: full match & not intersect ncRNAs:
  grep -wxFf $original_dir$prefix* $prefix.piRNA.fasta |grep -vwFf $dir2/bed/intersect_ncrna/all.ID.txt | sort -u | awk '{print ">read"NR"|FM\n"$1}'> $prefix.pi_reads.fasta
###                          changed 25/4 - add && NR!= 1 to #7.2
  #7.2: in loci of origin merged seq:
  awk '{if ($5 != "NR1" && NR!= 1) {print $9}}' $prefix.exp.bed | sort -u | grep -vwFf $prefix.pi_reads.fasta | awk '{print ">read"NR"|MR\n"$1}' >> $prefix.pi_reads.fasta
  #7.3: have 1-2 mm & 0 loci diff to seq:
  awk '{if ($5 == "NR1" && $8 == 0 ) {print $9}}' $prefix.exp.bed | sort -u | grep -vwFf $prefix.pi_reads.fasta | awk '{print ">read"NR"|MM\n"$1}' >> $prefix.pi_reads.fasta
  
  # 8: get seq pi fasta:
  grep -A 1 FM $prefix.pi_reads.fasta | grep -v ">" | awk '{if (NF!= 0){print ">piRNA_seq_1."NR"\n"$1}}'> $prefix.pi_seq.fasta
###                          changed 25/4 - add  &&  $1!="seq_f" 
  cut -f 4 $prefix.exp.bed | sort -u | grep -vwFf $prefix.pi_seq.fasta | awk '{if (NF!= 0 &&  $1!="seq_f" ){print ">piRNA_seq_2."NR"\n"$1}}' >> $prefix.pi_seq.fasta
  
  rm -rf $prefix.fitd.bed $prefix.unmtch2ori.txt $prefix.piRNA.bed $prefix.intrsct.bed tmp_$prefix.intrsct.txt 
  ##keep 4each infile, 5 files: $prefix.piRNA.fasta, $prefix.unmtch2ori.fasta , $prefix.exp.bed, $prefix.pi_reads.fasta, $prefix.pi_seq.fasta
  
  echo "done 1-loop" 
  date          
  
done

cat *pi_seq.fasta | grep -v ">" | sort -u | awk '{if (NF!= 0) {print ">piRNA_seq"NR"\n"$1}}' > piRNA_all_seq.fasta

# 9: get nonpiRNA reads
cat *pi_reads.fasta | grep -v ">" | sort -u | awk '{if (NF!= 0) {print ">piRNA_read"NR"\n"$1}}' > piRNA_all_reads.fasta
grep -vwxFf piRNA_all_reads.fasta $original_dir* | awk 'FS=":" {NR = NF ; {if (! /seq/) {print $NR}}}' |sort -u | awk '{if (NF!= 0){print ">nopiRNA"NR"\n"$1}}' > nopiRNA_all.fasta
# 10: seqmap
/gpfs0/alal/users/amitgefe/tools/seqmap-1.0.12-linux-64 2  piRNA_all_seq.fasta nopiRNA_all.fasta  seqmap_seq.out  /output_all_matches


#11: get_uncertin from seqmap_seq
/gpfs0/alal/projects/utils/python3.8 $tools_dir/get_uncertin.py
grep -vwxFf tmp_uncertin_pi.fasta piRNA_all_seq.fasta | grep -v seq_f | grep -v ">" -B 1 | awk '{if (/'--'/) {next}}{print $0}' > final_piRNA_seq.fasta 
###                          changed 1/5 - add awk '{if (length($0)< 5) {next}}{print  $0}'| grep -v '>\|N' -B 1  |
grep -vwxFf tmp_uncertin_nopi.fasta nopiRNA_all.fasta | awk '{if (length($0)< 5) {next}}{print  $0}'| grep -v '>\|N' -B 1  | awk '{if (/'--'/) {next}}{print $0}' > final_nopiRNA.fasta

#12. make fold files.
/gpfs0/alal/projects/utils/ViennaRNA-2.4.14/bin/RNAfold -i final_piRNA_seq.fasta  > final_piRNA.fold
find -name "*.ps" -delete
/gpfs0/alal/projects/utils/ViennaRNA-2.4.14/bin/RNAfold -i final_nopiRNA.fasta  > final_nopiRNA.fold
find -name "*.ps" -delete

#13. make csv file with labels
/gpfs0/alal/projects/utils/python3.8 $tools_dir/Vienna_2_kmers_dist.py -k 3 -i final_piRNA.fold -i final_nopiRNA.fold 

gzip $original_dir*
echo "done" 
date
