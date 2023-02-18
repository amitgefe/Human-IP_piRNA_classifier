#!/usr/bin/python
#filter sam file by seq length & min expression (count/ row[0][:row[0].find(':')+1])
import csv
import sys

def ReverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])

###                          changed 24/4 - remove +1 mini(11)&maxi(13) and -1 for seq.(29)
mini = int(sys.argv[1]) #chanegs to 22 #original 26

maxi = int(sys.argv[2]) #chanegs to 38 #original 36

min_expression = int(sys.argv[3]) #expression. pilfer set to 100, am_pipe set for 10, IP test set to 1

csvin = csv.reader(sys.stdin, delimiter="\t")

for row in csvin:
  if not row[0][0] == "@":
    f=row[0].split(":")
    row[0] = f[1]
    seq = row[9]  
    if len(row[9])>=mini and len(row[9])<=maxi and int(f[0])>= min_expression: #expression. pilfer set to 100, am_pipe set for 10, IP test set to 1
            strand = "+"
            if (int(row[1]) & (0x10)):
                seq = ReverseComplement(seq)
                strand = "-"
            #seq = seq[1:] #remove first char added by bowtie -- so the reads identical to the original fastq. 
            print (row[2] + "\t" + row[3] + "\t" + str(int(row[3])+len(seq)-1) + "\t" + seq + "::PU\t" + str(f[0])  + "\t" + strand) 
            
            
# example input line:
# 1:seq000016     0       chr15   62251025        255     21M     *       0       0       TGTGCTGTTTCTGCGGGGAGA   IIIIIIIIIIIIIIIIIIIII   XA:i:0  MD:Z:21       NM:i:0  XM:i:2

# example output line:
# GL000008.2      13130   13158   CATCCTCTCCAGCACCTGTTGTTTCCTGA::PU       1       +
