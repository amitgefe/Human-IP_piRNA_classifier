#!/usr/bin/python
from itertools import product


def get_kmer_dist(seq, k):
    '''
    The function get_kmer_dist get as an input a DNA representative string 
    with lower-case letter represents unpaired nucleotide and an upper-case represents
    paired nucelotide.
    The second input is an integer represents the length of the k-mers. 
    The function computes the relative abundance of each k-mer and return a distribution.
    
    Input:
    ------
        seq: string
        ----------
            A legitimate DNA representative string 
            
        k: int
        ------
            The desired length of k-mers 
    
    Output:
    -------
    List of floats represents the distribution of all possible k-mers in seq.
    '''
    try:
        # Get a dictionary of all possible kmers from 8 possible DNA bases (lower and upper case)
        kmers_count = {''.join(kmer): 0 for kmer in product(['a', 't', 'c', 'g', 'A', 'T', 'C', 'G'], repeat=k)}
        # Iterate over seq indices to generate kmer while updating kmer count in kmers_counts dictionary
        for i in range(len(seq)-k-1):
            kmers_count[seq[i:i+k]] += 1 
        # Sum up all kmers counts
        total = sum(kmers_count.values())
        # Generate a distribution of kmers instances in seq
        kmers_dist = [count/total for count in kmers_count.values()]
        return kmers_dist
    except Exception as err:
        raise  err(seq)



#x = get_kmer_dist('ctgtacAGCctcctaGCTttcc', 3)


#origin seq:
#CUGUACAGCCUCCUAGCUUUCC
#......(((......)))....         

#seq representative string:
#ctgtacAGCctcctaGCTttcc

