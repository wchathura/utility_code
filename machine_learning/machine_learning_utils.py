# This script contains helper functions to create test and traning datasets and to visualizing results
from collections import defaultdict
import pandas as pd
import numpy as np
import re

def extract_upstream(feature_file: pd.DataFrame, genome_dict: dict, up_down_len = [2000,0]):
    
    '''This function is used to extract specificed regions from up and down of the trnascription 
       Start site. 
       Inputs :
            "feature_file" is a cvs file where columns in the file are 
                                        0: chromosome id
                                        1: strand
                                        2: gene_id
                                        4: start_coordinate
                                        5: end_coordinate
            "geneome_dict" is a dictionary. Ecah key is a chromosomeID and value is a biopython "SeqRecord" object
                                        
                                        ex:genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

            "up_down_len" is a list specifying what coordinates to extract
                                        first element indicates upstream limit to extrct
                                        second element indicates downstrem limit to extract

                                        for ex: [2000,0] extract upstream 2000 base pairs from TSS
                                                [2000,150] extract upstream 2000 base pairs from TSS and downstrem 150 base pairs
        Output a dictionary. Each key in the output is a gene id. Value is a list with three elemnets. 0: chromosoeme, 1: strand, 2: sequence

    '''
    ## iterate through the feature file and extact chromosome id, strand, gene_id, start_coordinates and end_coordinates
    
    gene_dict=defaultdict(list)
    for f in feature_file.iterrows():
        # extarct chromosome id, strand, gene_id, start_coordinate, end_coordinate
        chro, strand, gene_id, start, end = f[1][0], f[1][1], f[1][2], f[1][4], f[1][5]
        
        ## convert the gene id to uppercase letters
        gene_id=gene_id.upper()
        # extract limits of up and downstrem regions 
        up_lim=up_down_len[0]
        down_lim=up_down_len[1]
        
        if strand=="+":
            #if gene is in the plus strand substract bases to get the upstrem region
            up_start=start+down_lim
            up_end=start-up_lim

            upstream_seq=genome_dict[chro][up_end:up_start]
            gene_dict[gene_id].append(chro)
            gene_dict[gene_id].append(strand)
            gene_dict[gene_id].append(str(upstream_seq.seq))

        if strand=="-":
 
            up_start=start-down_lim
            up_end=start+up_lim
            upstream_seq=genome_dict[chro][up_start:up_end]
            upstream_seq=upstream_seq.reverse_complement()

            gene_dict[gene_id].append(chro)
            gene_dict[gene_id].append(strand)
            gene_dict[gene_id].append(str(upstream_seq.seq))
    # in gene_dict key is a gene id, value is a list.  0 th element of the value list contains
    # chromosome id, 1 elemnet contains strand and the 2 element contains sequence    
    return gene_dict



def one_hot_encode(gene_seq: str):
    
    '''
        based on #https://stackoverflow.com/questions/34263772/how-to-generate-one-hot-encoding-for-dna-sequences
        Given a nucleotide sequence onehot encode it
        Input nucleotide sequence of arbitrary length
        Output one hot encoded 2D list 
    
    '''
    # define the encodings
    encoding_dict={0:[1,0,0,0,0],1:[0,1,0,0,0],2:[0,0,1,0,0],3:[0,0,0,1,0],4:[0,0,0,0,1]}
    
    gene_seq=gene_seq.upper()
    replaced_seq=gene_seq.replace("A",'0').replace("T",'1').replace("G",'2').replace("C",'3')
  
    # replace characters other than A, T, G, C as 4
    replaced_seq = re.sub(r"[A-Z]", "4", replaced_seq)
    int_s=[encoding_dict[int(s)] for s in replaced_seq]
    
    return int_s


def create_datasets(df: pd.DataFrame, seq_dict: dict, class_mapping: dict ,len_of_seq: int):
    '''
        This function is used to create the traning and test datasets for model traning 
        Inputs : df is a datafrme. This dataframe should have a "transcript_id" column which is the name of the 
                 transcript or the gene and "cat" column indicating class of that input for example in a binary classification 
                 case cat column could have 0 and 1 indicating the class of data point

                 seq_dict: Output a dictionary. Each key in the output is a gene id. Value is a list with three elemnets. 0: chromosoeme, 1: strand, 2: sequence
                           This dictionary is the output of the extract_upstream function

                 class_mapping: A dictionary where each key is the output class and a value is the encoding for that class
                                ex: {0:[0],1:[1]}
                                    {"cat": [0], "dog": [1]}
                 len_of_seq : lenght of the nuclotide sequence
        Outputs : One hot encoded array of nucleotide sequences and their targets
                
    '''
    
    X=[]
    Y=[]
    
    ## iterate throuh the df and extract the transcript id and the class it belongs to
    for g, ca in zip(df['transcript_id'],df['cat']):
        try: 
            ## upper is used here because in the dictionary all the geneIDs are in upper case
            gene_seq=seq_dict[g.upper()]
        except KeyError:
            "Gene id is not present in the input sequence dictionay"
            
        ## make sure correct output is there. 
        if len(gene_seq)==3:

            ## use the function one_hot_encode defined earlier to convert the string into one hot encoded vectors
            gene_encoded=one_hot_encode(gene_seq[2])
            if len(gene_encoded)==len_of_seq:
 
                X.append(gene_encoded)
                Y.append(class_mapping[ca])

    return np.array(X), np.array(Y)