import pytest
import numpy as np
import pandas as pd

from scsingscore.scsingscore import score
from scsingscore.scsingscore_ms import _ms_sing

DF = pd.DataFrame({
        'g1': [2, 1, 1, 0, 0],
        'g2': [1, 1, 0, 0, 0],
        'g3': [0, 0, 0, 1, 0],
        'g4': [0, 0, 0, 1, 0],
        }
    ).T


def test_singscore():
    up_gene = ['g1', 'g2']

    # note that we would usually remove the non-expressed genes.
    # but its fine for testing

    s = score(up_gene, DF, down_gene=False, norm_method='standard',
              norm_down=0, full_data=False, centering=True)

    expected = pd.DataFrame(
        {'total_score': {0: 0.375, 1: 0.25, 2: 0.125, 3: -0.25, 4: -0.25}}
    )
    assert np.all(expected == s)


def test_mssingscore():
    gset = ['g1', 'g2']

    # note that we would usually remove the non-expressed genes.
    # but its fine for testing
    expected = {0: 0.375, 1: 0.25, 2: 0.125, 3: -0.25, 4: -0.25}
    for i in range(len(DF.columns)):
         assert _ms_sing(gset, DF[i], norm_method='standard')['total_score'] == expected[i]
