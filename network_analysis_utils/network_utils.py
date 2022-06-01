import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import pearsonr, spearmanr
import pickle 

def openPickle(pk_file: str):
    """
        Open a pickle file
        Parameters:
                Location of the file 
        Returns:
              Opened pickle file    
    """
    with open(pk_file, "rb") as fi:
        return pickle.load(fi)


def writePickle(file_to_write: object, file_name: str):
    '''
        Write a object to a pickle file
        Parameters:
            file_to_write: File object that needs to be written 
            file_name: Name of the file to be writtern
    
    '''
    with open(file_name ,'wb') as handle:
        pickle.dump(file_to_write, handle, protocol=pickle.HIGHEST_PROTOCOL)


def calculteSigend_coor(coor_df: pd.DataFrame):
    """
        Calculate unsigned correlation 
        Based on "Weighted Network Analysis_ Applications in Genomics and Systems Biology-Springer-Verlag New York (2011)"
                  Page 99 
        
        Parameters: 
                Dataframe with correlation values. Values vary between -1 and +1
        Returns: 
             Dataframe with signed correlation values
             Correlation values will be in the range of 0 to 1
             -1 will be scaled to 0 and 1 one will be the same
    """
    
    # make sure dataframe is symmetrical
    assert coor_df.shape[0] == coor_df.shape[1]
    
    np.fill_diagonal(coor_df.values, 0)
    return (coor_df*0.5)+0.5


def calculate_coor_values(df: pd.DataFrame):

    '''
        Calculate the correlations between each columuns of a given dataframe
        This function uses the scipy.stats pearsonr, spearmanr functions to calculate
        the correlations

        Parameters :
                 samles X measurments datafrmae
                 ex: samples X expression of genes
        
        Returns : Two dataframes. First df contains correlations between each of the columns
                  Second datafrmae contains p-values for each correlation 
                  Both dataframes are of dimensions measurments X measuremnts 
    '''
    # code adopted from https://stackoverflow.com/questions/25571882/pandas-columns-correlation-with-statistical-significance
    df = df._get_numeric_data()

    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    coors = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1], 4)
            coors[r][c] = round(pearsonr(df[r], df[c])[1], 4)
    return coors,pvalues 


def getEdgeProperties(corr_df: pd.DataFrame, cut_off = 0.8):
    
    """ 
        Parameters: 
            2D numpy array of the correlation matrix
            dimensions N X N (for gene number)
            Iterate through each of the gene and calculate the connected genes and their strenghts
            
        Returns:
            2D dictionary:
            example : res_dict.keys() -> ["connections","strenghts]
            res_dict["connections"] is a dictionary where each key is a gene id and a 
            value is a list of nodes where the key is connected
            res_dict["connections"] -> {"gene1": ["gene2","gene3"],"gene2": ["gene7","gene19"]}
            res_dict["strenghts"] is a dictionary where each key is a gene id and a 
            value is a list of correlation strenght for each of the connected nodes
            res_dict["strenghts"] -> {"gene1": [0.8,0.85],"gene2": [0.92,0.85]}

    """
    con_net=defaultdict()
    strenght_net=defaultdict()
    
    ## set the diagnoal of the arry to -1
    np.fill_diagonal(corr_df.values, -1)
    

    res_dict={}
    i=0
    for idx_,f in zip(corr_df.index,corr_df.to_numpy()):
        
        con_net[idx_] = np.where(f >= cut_off)[0]
        strenght_net[idx_] = f[f >= cut_off]
    
        if i % 5000 == 0:
            print(i, "Of genes processed")
        i=i+1
    res_dict["connections"] = con_net
    res_dict["strenghts"] = strenght_net
    print("Done")
    return res_dict


def calculateNumberofConnections(net_dict: dict):
    ''''
        This function calculates the number of connections for each gene
        Input is a dictionary where key is the  gene id and the value is a list of nodes that are connected 
        to that key

        if the function getEdgeProperties has been called before we can pass the dictionary stored in the "connections"
        key

        for ex we can use the this function as below: 
        
            test_net = getEdgeProperties(edge_strenght_list: np.array, cut_off=0.8):
            num_of_connections=calculateNumberofConnections(test_net["connections"])

        Parameters: 
            net_dict : dictionary {"gene_id": [list of connected nodes]}

        Returns:

    '''
    num_of_connections=defaultdict()
    # iterate though the input connections dictionary
    # calculate the number of connected nodes as the lenght of the list
    for k,v in net_dict.items():
        num_of_connections[k] = len(v)
        #print(k,v)
    return num_of_connections

def calculateStrenghtOfconnections(net_dict: dict):
    
    '''
       Calculate the connection strength of a given node. Input is a dictionary where key is a gene id and 
       values are a list of correlation strengths for that perticular node
       For a given node sum up all the correlation strengths for a given node
       for ex: if for a node "gene1" strenghts are as follows, 0: 0.5 0.2 0.8 
       The function will return "gene1": 1.5 (sum up the connection strenghts)
        
        if the function getEdgeProperties has been called before we can pass the dictionary stored in the "strenghts"
        key

        for ex we can use the this function as below: 
        
            test_net = getEdgeProperties(edge_strenght_list: np.array, cut_off=0.8):
            network_strengths=calculateNumberofConnections(test_net["strenghts"])

    '''
    
    strenght_dict=defaultdict()
    for k,v in net_dict.items():
        strenght_dict[k] = sum(v)
        #print(sum(v))
    
    return strenght_dict


    
def convertToDeepwalk(adj: pd.DataFrame(),cutoff = 0.8):
    '''
       This is a helper function to convert adjacency matrix as input to the deepwalk algorithim (https://github.com/phanein/deepwalk)
       Adjacency matrix is convertd into adjacency list
       Output of this could be used to create the run the deepwalk algorithim with '--format adjlist' parameter

       Input is a adjacency matrix and the cutoff to use determine if two nodes are connected
       If the correlation strenght of the two nodes are higher or equal to the cut-off two nodes are considered 
       connected 
    '''
    i = 0
    deepWalk_dict=defaultdict()
    # iterate through the columns of the adj matrix, each column is a gene
    for f in adj.columns:
        # Extract the genes that have correlation more or eqal to the specified cut-off
        # use the gene as the key and append the set of gene as a list
        # Here we remove the current geneID since correlation of gene with itself is always 1

        deepWalk_dict[f] = set(adj[adj[f] >= cutoff].index.tolist())-set([f])
        i+=1
        if i%1000 == 0:
            print(i)
        
    ## return the dictionary as the output
    ## example output for "gene1" key is as follows
    ## deepWalk_dict["gene1"]=["gene2","gene3",gene5","gene7]
    return deepWalk_dict


def writeDeepWalkFormat(toWriteDict: dict, FileName: str):
    """
        This function takes the output dictionary from the function "convertToDeepwalk" and
        writes the adjacency list in to a text file. Inputs are output dictionary form convertToDeepwalk
        and a string specifying the output file  
        Output of this function could be used as an input to the DeepWalk algorithim
        Input is and adjecny 
    """
    ## open the file name
    with open(FileName,'w') as f:
        # iterate through each item in the dictionary 
        for k,v in toWriteDict.items():
            
            ## create a list by combining the key and the list 
            to_write = [k]+list(v)
            ## convert the list to a string joined by a spce
            LineToWrite = " ".join([str(kk) for kk in to_write])
            ## Write each line with a newline character at the end
            f.write(LineToWrite+"\n")


def formatToStellerGraph(conNet: dict) -> pd.DataFrame:
    
    '''
        This fucntion is used to create a pandas Df in the StellarGraph
        input format. 
        Input is dictionay. Input should have two keys "connections" and "strenghts"
        Each value for the "connections" and the "strenghts" is also a dictionary 
        These dictionary hold the information about the network 
        each key is a node and the value is a list. for connections key values are a list of nodes
        that are connectd to that particular node. for strenghts key, values are a list of connection strengths 
        for each node
    '''
    
    try: 
        connections = conNet["connections"]
        strenghts = conNet["strenghts"]
        
    except KeyError as e:
        print  ("'connections' or 'strenghts' key is not present in the dictionary")
    
    n1, n2, edge_strnght = [], [], []
    
    # iterate through each of the items in the net dictionary
    for i,(k, v) in enumerate(connections.items()):
        edge_values = strenghts[k]
        for vv,ed in zip(v,edge_values):
            #print(k, vv, ed)
            n1.append(k)
            n2.append(vv)
            edge_strnght.append(ed)
    
    graphDf = pd.DataFrame({"source":n1,"target":n2,"weight":edge_str})
    
    return graphDf 