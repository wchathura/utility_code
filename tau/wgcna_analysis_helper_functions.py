import pandas as pd
import numpy as np
import pyreadr
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union


def mergeEigWithSampleInfo(eig_expression: pd.DataFrame, sample_info: pd.DataFrame):
    """
    Merge eigen gene expression in replicates with the sample information table based
    on the sample name. Then get the mean of eigen gene expression in each or the replicate groups
    order the final output based on the order column and drop the order column
    """
    merged = pd.merge(
        eig_expression, sample_info, how="left", left_index=True, right_on="sample"
    )
    merged = merged.groupby("rep_group").mean()
    merged.sort_values(by="order", inplace=True)

    merged = merged.rename(columns=lambda x: x.split("ME")[-1])
    merged.drop("order", axis=1, inplace=True)

    return merged


def getExpressionOf_modules(
    expression_df: pd.DataFrame,
    mod_id: Union[int, str],
    module_id_column="netModuleLabels",
) -> pd.DataFrame:
    """
    Given the module ID return a datafrmae contaning the expression of all the genes in that modules

    Inputs: expression_df dataframe is [genes X sample] dataframe (from wgcna output)
            This dataframe should have a column to idedntfiy which column indicates the
            gene module relationship. This is set by the parameter "module_id_column
            Default is "netModuleLabels" which is WGCNA defaut identifier
            Mod id could be int or str depending on which moudle naming standard is follwed

     Output: return a dataframe of expression where
    """

    module_genes = expression_df[expression_df[module_id_column] == mod_id]
    return module_genes.reset_index()


def getEigenGeneExpression(eigen_expression_df: pd.DataFrame, mod_id) -> pd.DataFrame:
    """
    Inputs : eigen gene expression dataframe of ["samples" X eigenGene]
    Output : return a series of selected eigne gene expression
    """
    return eigen_expression_df[mod_id]


def plotWgcnaEXp(
    expDf: pd.DataFrame,
    eigDf: pd.Series,
    plt_size=(10, 5),
    method=0.95,
    plot_type="bar",
):
    """
    Plot eigne gene expression in samples and the expression of genes in the samples in
    the same plot using indpendnet Y-axis. (two y-axes have two seperate scales)
    Inputs :
            expDf: Dataframe contaning expression data of genes [Genes X samples]
                   This dataframe need to have a column named index that identifies the gene id
                   Also this dataframe should have a column named "netModuleLabels" whcih indicates
                   the gene module id of that particular gene

            eigDf:  Series contaning eigen gene expression of samples.

            method:  This paramters indicates wihch error bars to use default is to use 95 confidance interval

            plot_type: This parameter indicates which type is plot to plot: options could be line and bar
    """
    try:
        ## data needs to be in tje long format for sns.bar function
        melted_exp = expDf.melt(id_vars=["netModuleLabels", "netModuleColors", "index"])

    except KeyError as e:
        pass

        try:
            melted_exp = expDf.melt(id_vars=["netModuleLabels", "index"])

        except KeyError as e:
            print(e)

    # print(melted_exp)
    sns.set(font_scale=1.2)
    sns.set_style("ticks")

    fig, ax1 = plt.subplots(figsize=plt_size)
    color = "tab:red"
    ax1.set_xlabel("Samples")
    ax1.set_ylabel("eigneGenes", color="salmon")
    ax1.bar(np.arange(len(eigDf)), eigDf.values, color="salmon")
    # ax1.plot(np.arange(len(eigDf)),eigDf.values,color="salmon")
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = "teal"
    ax2.set_ylabel("Expression", color=color)  # we already handled the x-label with ax1

    if plot_type == "bar":
        ax2 = sns.barplot(
            data=melted_exp,
            x=melted_exp["variable"],
            y=melted_exp.value,
            ci=95,
            alpha=0.5,
            color="teal",
        )

    if plot_type == "line":
        ax2 = sns.lineplot(
            data=melted_exp,
            x=melted_exp["variable"],
            y=melted_exp.value,
            ci=95,
            alpha=0.5,
            color="teal",
        )
    ax2.tick_params(axis="y", labelcolor="teal")

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.xticks(np.arange(len(eigDf.index)), eigDf.index)
    ax1.tick_params(labelrotation=90)
    # plt.xticks(rotation=90)
