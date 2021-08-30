from anndata import AnnData
import numpy as np
from scsingscore.smoothing import get_smoothing_matrix
from scipy import sparse

def test_get_smoothing_matrix():
    """
    assert the the rows of the smoothing matrix are normalized to 1
    """
    ncells = 4
    ngenes = 3

    adjacency = np.array([
        [0, 0.1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]])

    adata = AnnData(np.random.rand(ncells,ngenes))
    adata.obsp['distances'] = sparse.csr_matrix(adjacency)


    S = get_smoothing_matrix(adata, mode='adjacency')
    np.testing.assert_allclose(S.A.sum(1), 1)


    adata.obsp['connectivities'] = sparse.csr_matrix(adjacency)
    S = get_smoothing_matrix(adata, mode='connectivity')
    np.testing.assert_allclose(S.A.sum(1), 1)


    """
    make sure it throws an error when the matrix has diagonal elements
    """
