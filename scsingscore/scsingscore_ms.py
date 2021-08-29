"""
MS version of scsingscore.py
"""
import numpy as np
from scipy import sparse
import pandas as pd
import scsingscore.scsingscore as si
import tqdm
import statsmodels.robust.scale
import scanpy as sc
from scsingscore.smoothing import nn_smoothing


def sc_score(
        adata,  # anndata containing single cell rna-seq data
        noise_trials,  # number of noisy samples to create, integer
        samp_neighbors,   # number of neighbors to sample
        gene_set,  # the gene set of interest
        mode='average',  # average or theoretical normalization of scores
        ):

    # NOTE: this is cells x genes
    smoothed_matrix = nn_smoothing(adata.X, adata, 'connectivity', samp_neighbors)
    # for easier handling with gene names
    smoothed_adata = sc.AnnData(smoothed_matrix, obs=adata.obs, var=adata.var)
    """
    since we're doing all cells at the same time now,
    the following gets probelmatic (the df kills the sparsity)
    Ideas:
    - loop over cells, do the scoring
    - batch the cells, i.e. create a df with 100 cells (puling them into mem) and score those in one go
    """
    # # then we subset it to only the genes with counts
    # df = pd.DataFrame(gene_mat.T.A, index=adata.var.index, columns=adata.obs.index)
    # gdx = df.sum(1) == 0
    # df = df.iloc[gdx,:]
    # df.columns = ['gene_counts']
    all_scores = _score_one_by_one(gene_set, smoothed_adata, noise_trials, mode='average')
    return all_scores


def _score_one_by_one(gene_set, smoothed_adata, noise_trials, mode='average', ):
    scores_across_cells = []
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):
        gene_mat = smoothed_adata.X[cell_ix]
        # then we subset it to only the genes with counts
        _, gdx, _ = sparse.find(gene_mat)
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=smoothed_adata.var.index[gdx])
        df.columns = ['gene_counts']

        # FROM HERE: parallel computing of scores across gene sets.
        # or parallelize on each column of dataframe

        if mode == 'average' and noise_trials > 0:
            # add some noise to gene counts.. create a n numbers of examples
            df_noise = si.add_noise(df, noise_trials, 0.01, 0.99)  # slow part .. fixed
        else:
            df_noise = df
        # score the neighborhoods
        s = si.score(up_gene=gene_set, sample=df_noise, norm_method='standard', full_data=False) # standard workin gbetter here than theoretical

        avg_score = s.mean()
        scores_across_cells.append(avg_score)

    return np.array(scores_across_cells).flatten()
