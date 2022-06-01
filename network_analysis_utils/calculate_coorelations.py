import argparse
import pandas as pd
import numpy as np

from scipy.stats import pearsonr
import pandas as pd
import network_utils as net_utils

## calculate coorelations

group_dict = { "tre_groups" :[ '4DST', 'SpSD12T', 'SpSD8T','SpRD8T', 'SpRD12T','SpRMu3HT','SpRMuD1T','SpRMuD3T','SpSMu3HT',
 'SpSMuD1T', 'SpSMuD3T','SpL567W5T','YLT','SpST7T','SpFW6T','SpFW8T','SpW2SQT','SpW6SQT','SpSeeT'],

	      "con_groups" : ['2DSC','4DSC','SpSD8C','SpRD8C','SpSD12C','SpRD12C','YLC','SpL567W5C','SpST7C','SpFWC','SQC',
			'SpSeeC','SpRMuC','SpSMuC'] 
              }

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Coorelation calculation')
    parser.add_argument('-f', type=str, required=True, default=False)
    parser.add_argument('-o', type=str, required=True, default=False)
    parser.add_argument('-r', type=str, required=False, default="pearson")
    args = parser.parse_args()
    

    multi_df=pd.read_csv(args.f,index_col=[0,1])
    tpm_mean = multi_df.groupby("rep_group").mean()
    
    for name, gr in group_dict.items():
        df = tpm_mean[tpm_mean.index.isin(gr)]
        cor, p_val = net_utils.calculate_coor_values(df)
        print(args)

        cor.to_csv(args.o + "_"+ name +"_cor_out.csv")
        p_val.to_csv(args.o + "_"+ name +"_p_out.csv")

