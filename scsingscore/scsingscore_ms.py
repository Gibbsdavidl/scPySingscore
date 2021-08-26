"""
MS version of scsingscore.py
"""
import numpy as np
from scipy import sparse
import pandas as pd
NN_DISTANCE_KEY = 'distances'

# TODO test: should always sum to 1
# multiplying should leave a "one-vector" still sum to one


def get_smoothing_matrix(adata, mode):
    if mode == 'adjacency':
        # TODO: this does not contain the cell itslef?!
        # TODO: make sure the neigbhours of cell i are in A[i] (not A[:,i])
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

    return sparse.csr_matrix((newvals, (newrows, newcols)))  # TODO choose matrix format as X


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
    gene_mat = nn_smoothing(adata.X, adata, 'connectivity', samp_neighbors)

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
