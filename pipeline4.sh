#!/bin/bash
##      preprocessing of raw data (fasta/q files)      ##
## from collapse within files to genome alighment, length filtretion, remove known ncRNA and overlap genes.##
#update of /gpfs0/alal/users/amitgefe/am_pipe/supp/pipeline1.sh & /gpfs0/alal/users/amitgefe/IP_data/preprocess/get_piRNA_can.sh
# origin from pipeline1.sh, Collapse to  "filter sam" (length) include
# origin from get_piRNA_can.sh,  ncrna intersect to putative/gene_subtract/

#use qsubs_split.sh

tools_dir="/gpfs0/alal/users/amitgefe/am_pipe/supp/"
genome_dir="/gpfs0/alal/users/amitgefe/genome_chm13v2/"

###                          changed for train set                          ###
out_dir="/gpfs0/alal/users/amitgefe/preprocess_chm13v2.0/" #all 266 samples with new genome
##"/gpfs0/alal/users/amitgefe/IP_classifier/train/"  # for train (IP data)  OR /train_chm13v2.0/" for new genome
##"/gpfs0/alal/users/amitgefe/preprocess/" #for test (all data 266) 
##"/gpfs0/alal/users/amitgefe/preprocess_chm13v2.0/" #all 266 samples with new genome

###                          changed for train set                          ###
sample_dir="/gpfs0/alal/users/amitgefe/availble_data/" # for test (all data 266)
#"/gpfs0/alal/users/amitgefe/IP_data/test/" # for train (IP data) - *.amtiTrim.fastq
#"/gpfs0/alal/users/amitgefe/availble_data/" # for test (all data 266)

echo "sample_dir = ${sample_dir}"
echo "out_dir = ${out_dir}"

#for each sample
echo $prefix
echo $prefix >> $out_dir/prefix_list
date

echo "Running Collapse"
gzip -d -c  $sample_dir$prefix".fastq.gz" | /gpfs0/alal/projects/utils/python3.8 $tools_dir"pilfer_collapse.py" -i - -o $out_dir"collapse/"$prefix".fa" --format fastq

echo "Mapping to hg38"
gzip -d -c $out_dir"collapse/"$prefix".fa.gz" | /gpfs0/alal/projects/utils/bowtie-1.2.3-linux-x86_64/bowtie -p 8 /gpfs0/alal/projects/utils/bowtie-1.2.3-linux-x86_64/indexes/T2T-CHM13v2.0.genome -f - -a --best --strata -l 22 -n 1 -S $out_dir/hg38_sam/$prefix.sam 2> $out_dir/logs/bowtie1/"$prefix"_bowtie1.log

echo "filter sam"
###                          changed for train set                          ###
# filter for len(22-36)& min_expression >= 3. for train set  min_expression >= 1
# startswith int=  expression in origin fastq. in pilfer bed format. F 4 == without read unmapped 
/gpfs0/alal/projects/utils/samtools/samtools-1.10/samtools view -h -F 4 $out_dir"hg38_sam/"$prefix.sam | /gpfs0/alal/projects/utils/python3.8 $tools_dir"filter_sam.py" 22 36 3| sort -V -k1,1 -k2,2 > $out_dir"bed/filtered/"$prefix.bed
#to save memory:
/gpfs0/alal/projects/utils/samtools/samtools-1.10/samtools view --no-PG -bh -S -F 4 $out_dir"hg38_sam/"$prefix.sam > $out_dir"hg38_sam/"$prefix.bam 
rm $out_dir"hg38_sam/"$prefix.sam
  
echo "ncrna intersect"
#first output all loci contain any fregment of sncRNA loci, then get all the loci of the reads without any overlap.
#get all the reads intersect with known sncRNA out of putative bed file. by seq (not loci)
/gpfs0/alal/projects/utils/bedtools2/bin/bedtools intersect -nonamecheck -b $genome_dir"RNACentral.20.T2T.sncRNA.bed"  -a $out_dir"bed/filtered/"$prefix.bed -s -f 0.5 -wa > $out_dir"bed/intersect_ncrna/"$prefix.bed
  
cut -f 4 $out_dir"bed/intersect_ncrna/"$prefix.bed | sort -u  > $out_dir"bed/intersect_ncrna/"$prefix.ID.txt
  
echo "ncrna subtract"
grep -wvFf $out_dir"bed/intersect_ncrna/"$prefix.ID.txt $out_dir"bed/filtered/"$prefix.bed > $out_dir"bed/ncrna_subtract/"$prefix.bed
rm -rf $out_dir"bed/intersect_ncrna/"$prefix.bed
  
#merge by gemonic location. ($4-$6: total expression- XP, amount of different reads-NR, strand)
/gpfs0/alal/projects/utils/bedtools2/bin/bedtools merge -i $out_dir"bed/ncrna_subtract/"$prefix.bed -d -1 -s -c 4,5,6 -o count_distinct,sum,distinct |awk ' {$4 = "XP:" $5 "::NR:" $4 ; print; } '  | tr " " "\t" > $out_dir"putative/"$prefix"_merged.bed"
  
#genes subtract by location (not read)
echo "genes subtract"
/gpfs0/alal/projects/utils/bedtools2/bin/bedtools subtract -nonamecheck -b $genome_dir"CHM13.v2.0.gff3" -a $out_dir"putative/"$prefix"_merged.bed" -s -f 0.5 -A > $out_dir"putative/gene_subtract/"$prefix.bed

###                          changed 25/4 - skip & >=37)
# get sequences & filter length 
#awk 'OFS="\t" {$2=$2-1; if ($3-$2+1>=37) {next}} {print $0}' $out_dir"putative/gene_subtract/"$prefix.bed | /gpfs0/alal/projects/utils/bedtools2/bin/bedtools getfasta  -s -name+  -fi $genome_dir"chm13v2.0.genome.fa" -bed -  -fo $out_dir"putative/gene_subtract/"$prefix.fasta

echo "done"
date

