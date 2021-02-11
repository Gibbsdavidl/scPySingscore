
# one cell, one score

import scsingscore.scsingscore as si
import scanpy as sc

# get our gene sets
gd = si.read_gene_sets('F:\Work\cruk\data\epithelial_proliferation_geneset.gmt')
gene_set = gd['GO_EPITHELIAL_CELL_PROLIFERATION']  # 'HALLMARK_TNFA_SIGNALING_VIA_NFKB']

# read in more gene sets
gd = si.read_gene_sets('F:\Work\cruk\data\h.all.v7.2.symbols.gmt')
gs1 = gd['HALLMARK_TNFA_SIGNALING_VIA_NFKB']
gs2 = gd['HALLMARK_E2F_TARGETS']
#print(len(gd['HALLMARK_TNFA_SIGNALING_VIA_NFKB']))
print(gs1[0:5])
gs_list = [gs1, gs2]


q = sc.read_h5ad("F:\Work\cruk\data\\test_data.h5ad")

# what cell we'll look at
celli = 55
noise_trials = 3  # number of noise cells to make
num_neighbors = 33
samp_neighbors = 27
mode = 'average'

# 33% of time was on running neighbors
sc.pp.neighbors(q, n_neighbors = num_neighbors)

for celli in range(1,99):

    # 33% of time is on scoring call... can't change that much
    # 12% of time is to_list()
    onescore = si.one_score(q, celli, 
        noise_trials, num_neighbors, samp_neighbors, gene_set, 
        mode, False)
    print('done ' + str(celli))
    #print(onescore.total_score)


