import pandas as pd
import numpy as np
from sklearn.preprocessing import scale, normalize
from sklearn.decomposition import PCA 
import seaborn as sns
import umap


def run_pca(norm_df,columns=[],pca_compo=2):

    pca = PCA(n_components=pca_compo)
    pr_Compo= pca.fit_transform(norm_df.T)
    return pr_Compo

def run_umap(norm_df,columns=[],umap_compo=2):
 
    fit_umap = umap.UMAP(n_components=umap_compo)
    umap_fitted = fit_umap.fit_transform(norm_df)
    return umap_fitted

def convertCountsToTPM(counts: pd.DataFrame, gene_lens: pd.DataFrame) -> pd.DataFrame:
   
    """
        This function is used to normilzae gene expression counts to TMP value 
        # based on https://support.bioconductor.org/p/91218/
        Inputs: 
                counts dataframe : genes X samples
                effectiive length dataframe : genes X samples
        Output:
                Dataframe contaning TMP normilized counts 
    """

def calculateEffect_size(expDf: pd.DataFrame, effects_list: list)-> pd.DataFrame:
    
    '''
        Inputs: 
            Calulate the effect size given a expression dataframe and set of effects.
            Effect is defined as the change form contol to treatment.
            Treatment to Control ratio is caluclated ex treS1/ conS1
            
            expDf : pandas dataframe conaning the mean expression for samples. Columns are samples
                    rows are genes. 
            effects list : 2D list. each elemnt of the list contains two items which are used to 
                           calculate the effect. control should be the first element 
                           for ex  test_list = [[conS1, treS1],[conS2, treS2]]
    '''
    
    # Declare a new  pandas datafrmae
    effect_df = pd.DataFrame()
    try:
        # Check if the key values are present in the expression dataframe
        for comp_group in effects_list:
            #print(comp_group)
            # Control condition needs to be in the 0th index
            effect = expDf[comp_group[1]] / expDf[comp_group[0]]
            effect_df["__".join(comp_group)] = effect
    
    except KeyError as k:
        print(f"Keys {k} not found in the Expression Dataframe")
    
    return effect_df
