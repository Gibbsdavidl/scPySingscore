
import scanpy as sc
import scsingscore.scsingscore as si
from scsingscore.scsingscore_ms import sc_score_ms

gd = si.read_gene_sets('F:\Work\cruk\data\h.all.v7.2.symbols.gmt')
gs1 = gd['HALLMARK_TNFA_SIGNALING_VIA_NFKB']
gs2 = gd['HALLMARK_E2F_TARGETS']

q = sc.read_h5ad("F:\Work\cruk\data\\test_data.h5ad")
sc.pp.neighbors(q, n_neighbors=32)

# note that we would usually remove the non-expressed genes.
# but its fine for testing
expected = {0: 0.375, 1: 0.25, 2: 0.125, 3: -0.25, 4: -0.25}

x = sc_score_ms(
        adata=q,
        noise_trials=0,
        samp_neighbors=0,
        gene_set_up=False,
        gene_set_down=gs2,
        mode='average'
        )
print(x)