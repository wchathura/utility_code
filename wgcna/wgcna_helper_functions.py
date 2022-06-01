## helper functions for wgcna analyisis

import pandas as pd


def prepFor_Bingo(wgcna_df: pd.DataFrame,file_path: str):

    '''
        Prepare WGCNA output for bingo batch run
    '''
    batchBingo= open(file_path, 'w',encoding='utf_8')
    for gr in wgcna_df.groupby("netModuleLabels"):
        #print(gr[0],gr[1]["at_id"].tolist())
        batchBingo.write("%s\n" % (str(gr[0])+"_wgcna"))
        print(gr[0],len(gr[1]["gene_id"].tolist()))
        for gene in gr[1]["gene_id"].tolist():
            batchBingo.write("%s\n" % gene)
        batchBingo.write("%s\n" % "batch")

    batchBingo.close()