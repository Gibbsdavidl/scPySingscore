"""
MS version of scsingscore.py
"""
import numpy as np
from scipy import sparse
import pandas as pd
NN_DISTANCE_KEY = 'distances'
import scsingscore.scsingscore  as si
import tqdm


# TODO test: should always sum to 1
# multiplying should leave a "one-vector" still sum to one


def get_smoothing_matrix(adata, mode):
    if mode == 'adjacency':
        A = (adata.obsp[NN_DISTANCE_KEY] > 0).astype(int)

        # add the diagnoal, ie. the datapoint itself should be represented in the smoothing!
        A = A + sparse.diags(np.ones(A.shape[0]))

        # normalize to sum=1 per  row
        row_scaler = 1 / A.sum(axis=1).A.flatten()
        normA = A @ sparse.diags(row_scaler)
        return normA

    # actually works exactly the same; could just switch out the obsp key
    elif mode == 'connectivity':
        A = (adata.obsp['connectivities'] > 0).astype(int)
        # add the diagnoal, ie. the datapoint itself should be represented in the smoothing! Note that the max connectivity == 1
        A = A + sparse.diags(np.ones(A.shape[0]))
        # normalize to sum=1 per  row
        row_scaler = 1 / A.sum(axis=1).A.flatten()
        normA = A @ sparse.diags(row_scaler)
        return normA
    else:
        raise ValueError(f'unknown mode {mode}')


def random_mask_a_nn_matrix(X, nn_to_keep):
    """
    for a nearest neighbour matrix (i.e. each row has N entries)
    subsample the neighours (setting some entries per row to 0)

    this is pretty slow, maybe there's a better way...
    """
    newrows = []
    newcols = []
    newvals = []
    for i in range(X.shape[0]):
        _, col_ix, vals = sparse.find(X[i])
        rand_ix = np.random.choice(len(col_ix), size=nn_to_keep, replace=False)
        newrows.extend([i]*nn_to_keep)
        newcols.extend(col_ix[rand_ix])
        newvals.extend(vals[rand_ix])

    # TODO choose matrix format as X
    return sparse.csr_matrix((newvals, (newrows, newcols)))


def nn_smoothing(X, adata, mode, samp_neighbors):
    """
    :param X: data to smooth. (cell x feature) matrix
    :param adata: sc.AnnData, containing the neghbourhood graph
    :param samp_neighbors: consider all neighbours or a subsample of size samp_neighbors
    """
    assert X.shape[0] == adata.shape[0], "cell number mismatch"
    smoothing_mat = get_smoothing_matrix(adata, mode)

    if samp_neighbors > 0:
        # randomly set some edges/neighbours to zero
        # experimental and pretty slow!
        smoothing_mat = random_mask_a_nn_matrix(smoothing_mat, samp_neighbors)

    smooth_signals = smoothing_mat @ X

    return smooth_signals


def sc_score(
        adata,  # anndata containing single cell rna-seq data
        noise_trials,  # number of noisy samples to create, integer
        samp_neighbors,   # number of neighbors to sample
        gene_set,  # the gene set of interest
        mode='average',  # average or theoretical normalization of scores
        ):

    # NOTE: this is cells x genes
    full_matrix = nn_smoothing(adata.X, adata, 'connectivity', samp_neighbors)

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
    scores_across_cells = []
    for cell_ix in tqdm.trange(adata.shape[0]):
        gene_mat = full_matrix[cell_ix]
        # then we subset it to only the genes with counts
        _, gdx, _ = sparse.find(gene_mat)
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=adata.var.index[gdx])
        df.columns = ['gene_counts']

        # FROM HERE: parallel computing of scores across gene sets.
        # or parallelize on each column of dataframe

        if mode == 'average' and noise_trials > 0:
            # add some noise to gene counts.. create a n numbers of examples
            df_noise = si.add_noise(df, noise_trials, 0.01, 0.99) ## slow part .. fixed
        else:
            df_noise = df
        # score the neighborhoods
        s = si.score(up_gene=gene_set, sample=df_noise, norm_method='standard', full_data=False) # standard workin gbetter here than theoretical

        avg_score = s.mean()
        scores_across_cells.append(avg_score)
