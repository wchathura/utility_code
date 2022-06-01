import pandas as pd
import numpy as np
import sys
import argparse

def readDf(fileName: str):
    df=pd.read_csv(fileName)
    return df

def megreAdj_with_ModuleInfo(adj: pd.DataFrame(),modInfo: pd.DataFrame()):
    """ Merge the adjency matrix with the module info dataframe
        Module info dataframe has information about which gene belongs to which co-expressed cluster
    """
    # make sure modInfo df has moudle lable assingned
    try:
        dfSel= modInfo['bwnetModuleLabels','bwnetModuleColors']     
    except KeyError:
        print("Columns [bwnetModuleLabels','bwnetModuleColors] are not in" )
    
    ## merge adjency df and module names df
    merged=pd.merge(adj,dfSel,how="left",left_index=True,right_index=True)
    return merged

def extractGenesFromModules(modInfo: pd.DataFrame()):

    # make sure modInfo df has moudle lable assingned
    try:
        ## group genes by modules and return as a dict
        modGrouped=modInfo.groupby("bwnetModuleLabels").groups
        return modGrouped
    except KeyError:
        print("bwnetModuleLabels is not in the datafrme, Module lables are not assigned" )


def extractAdj_matrixForModules(adj: pd.DataFrame,moddict: dict())





def ()



def findTopHubGenes(modAdjDf: pd.DataFrame(),hub_num=15):
    """ From the adjancy matrix return the genes with the highest mean adjencey"""
    
    #check if matrix is a squre one
    assert modAdjDf.shape[0]==modAdjDf.shape[1]
    
    ## replace the diagnoal with 0, since diagnoal has always 1
    modAdjDf.values[[np.arange(modAdjDf.shape[0])]*2] = 0
    topHub=(modAdjDf.sum(axis=1))/(modAdjDf.shape[0]-1)
    #print(topHub.sort_values(ascending=False))
    return topHub.sort_values(ascending=False).head(n=hub_num)



def main():


if __name__ == "__main__":

    run(main)
    parser = argparse.ArgumentParser(description='Hub genes identfication')
    parser.add_argument('-f', type=str, required=True, default=False)
    parser.add_argument('-o', type=str, required=True, default=False)
    args = parser.parse_args()