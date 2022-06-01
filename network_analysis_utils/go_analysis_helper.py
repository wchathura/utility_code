import os
import glob
import pandas as pd

def format_bingo(dir_name: str,comparasion ="comp"):
    '''
        Given a set of bingo GO enrichment output file directory, extract the
        GO result and concatanate to produce a pandas dataframe

        Parameters :
            dir_name : string pointing to the bingo output files

        Returns : formated pandas dataframe, where all the bingo results are concatanated 
                  addtional column is added to indicate where the samples came from 

    '''
   

    # list the files in the directory 
    # files should be outputs from bingo with .bgo extension 
    # use the default bingo output format (.bgo extension) 
    # iterate through each of the bingo output files 
    print("jjjjj")
    #k=glob.glob(dir_name+"/*.bgo")
    #print(k)

    df_list = []
    for f_name in glob.glob(dir_name+"/*.bgo"):
        #f_name = "/*.bgo"
        i = 0 
        #print(f_name)
        with open (f_name) as go_res:
            for line in go_res.readlines():
                    i = i+1
                    #
                    # Here we look for a line which starts with Go-ID
                    # all the lines following this is formated as a pd.Dataframe
                    # use pd.read_csv function with skiprows to extrcat only the GO enriched results 
                    # --- may be find a another more efficient way later !! ---
                    if line.startswith("GO-ID"):
                        print(i, print(f_name.split("/")[-1]))
                        df = pd.read_csv(f_name,skiprows=i-1,sep="\t")
                        df[comparasion] = f_name.split("/")[-1]
                            
                        df_list.append(df)
   
    concat_df = pd.concat(df_list)
    
    return concat_df