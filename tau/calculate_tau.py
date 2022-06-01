import pandas as pd
import numpy as np
from scipy import stats


#def createMulti_index_df(tpmDf, colour_mapping_df):
#    """
#    Create create a multi index dataframe from a dataframe that has tpm values
#    and a datafrmae that has mapping from samples to rep groups
#    This is useful in calculating means, std etc
#    """
#    pass


def removeZeroExpGenes(df: pd.DataFrame) -> pd.DataFrame:
    """
    remove genes with zero expression in all the samples
    Dataframe should be Samples by genes matrix
    """

    exp_fil = df.drop(df.columns[df.sum(axis=0) == 0], axis=1)
    return exp_fil


def removeLowExpressedGenes(df: pd.DataFrame, expThresh=2) -> pd.DataFrame:
    """Group genes based on the replicate groups and calculate the mean
    Filter the genes where mean is less than 2 tpm
    """
    # return df[df.columns[(df.groupby(level=1).mean()>=2).any()]]
    sel_cols = df.columns[(df.groupby(level=1).mean() >= expThresh).any() == True]

    return df[sel_cols]


def calculate_tau(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the tissue specificty index according to (Yanai, Itai, et al., (2005)
    Values need to be destributed between
    Steps: Normalize the dataframe by maximum value.
    Steps :
            1) Divide the whole dataframe by maxmium value of the datafame.
            2) Substract 1 from all the values
            3) Sum across rows and divided by the number of (samples-1)
    """
    df_max_normlaized = df.apply(lambda x: x / x.max(), axis=0)
    df_one_substratced = df_max_normlaized.apply(lambda x: 1 - x, axis=0)
    tau_df = df_one_substratced.sum(axis=0) / (df_one_substratced.shape[0] - 1)
    return tau_df


def removeVariableGenes(
    ExpDf: pd.DataFrame, rep_number=3, z_score_thresh=2
) -> pd.DataFrame:

    """
    This function will filter the genes out the genes that are highly variable across biological replicates using the z-score of
    gene expression acoross replicates
    In more than half of all the samples in the set need to have z-score less than 2 in all replicates
    Input dataframe has to be a multi-index dataframe
    Second index should be the replicate goup
    first index is for the sample and the second index is for the replicate group

    """

    ## make sure each group has at least three replicates
    assert (ExpDf.groupby(level=1).count().iloc[:, 0] < 3).sum() == 0
    samples = np.ceil(ExpDf.groupby("rep_group").count().shape[0] / 2)

    # Calculate the mean of gene expression in each replicate group
    multi_mean = ExpDf.groupby(level=1).transform("mean")

    # Calculate the STD of gene expression in each replicate
    multi_std = ExpDf.groupby(level=1).transform("std")

    # Substrat mean from the expression dataframe
    mean_sub = ExpDf - multi_mean

    # Calculate the absolute Z-Score of each gene
    abs_z = (mean_sub / multi_std).apply("abs")

    # get the count of samples where the Z score of gene expression in tpm is less than or equal to 2
    # in at least in three replicates
    z_lessThan2 = abs_z.groupby(level=1).apply(
        (lambda x: (x <= z_score_thresh).sum() >= rep_number)
    )

    # filter out genes where the expresison of the gene is  variable in more than half of the samples

    filtered = ExpDf[z_lessThan2.columns[(z_lessThan2.sum() > samples)]]

    return filtered


def identify_significant_modules(
    eigen_df : pd.DataFrame, comp_list: list, p_value_cut_off=0.001
):
    """
    This function compares the mean of given set of conditions to
    the other set of conditions. ex: mean of ["x1","x2","x3"] to mean of ["y1","y2","y3","y4"]
    means are compared using individual two independent samples t-test for and setting
    equal_var=False
    """

    # select the eigen values corresponding to the input conditions
    try:
        sel_con = eigen_df.loc[comp_list]

    except KeyError as e:
        "Given samples are not present in the eigen gene expression dataframe"

    # select the eigen values corresponding to the input conditions
    oth_con = eigenDf[~(eigenDf.index.isin(comp_list))]

    # define an empty list to store the eigen geneIDs and p-values
    sig_list = []

    # go through the each column of the selected comprasions dataframe
    # coloum corrosponds to eigen gene cofficient for each comparasion
    # extract the corrosponding eigen gene from the dataframe to compare and perform the test
    for s in sel_con:
        test_res = stats.ttest_ind(sel_con[s], oth_con[s], equal_var=False)
        if test_res[1] <= p_value_cut_off:
            sig_list.append((s, test_res[1]))

    df = pd.DataFrame.from_records(sig_list, columns=["eigenID", "p_value"])

    return df
