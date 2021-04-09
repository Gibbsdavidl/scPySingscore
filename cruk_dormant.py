import pandas as pd
from scsingscore.scsingscore import *

def convertMouseSet(settext, hmd, seti):
    mouse_genes = settext[seti].split('\t')[1].split(',')
    human_genes = []
    for mi in mouse_genes:
        res0 = hmd.loc[hmd.MouseSymbol == mi.strip(), 'HumanSymbol'].to_list()
        if len(res0) > 1:
            print("WARNING: mouse gene has multiple matches: " + str(res0))
            human_genes += res0
        else:
            human_genes.append(res0[0])
    return(human_genes)


# dormancy gene sets from Julio et al.
settext = open('F:\Work\cruk\dormancy\dormancy_mouse_genelists.txt', 'r').read().split('\n')
# HMD gene orthologs from Jax 4/8/21
hmd = pd.read_csv('F:\Work\cruk\dormancy\HMD_HumanPhenotype.rpt', sep='\t',  index_col=False )
# cluster 9 markers
test_set = ['TPSAB1', 'TPSB2', 'CPA3', 'HPGDS', 'MS4A2', 'KIT', 'IL1RL1', 'LTC4S', 'SRGN', 'ALOX5AP', 'GATA2', 'FCER1G',
            'CD69', 'RGS1', 'VIM']

gene_set = convertMouseSet(settext, hmd, 0)

set_list = [test_set, gene_set]

i = 9
thisq = sc.read_h5ad('F:\Work\cruk\dormancy\data\clusters\esoatlas_c' + str(i) + '.h5ad')

# what cell we'll look at
noise_trials = 7  # number of noise cells to make
num_neighbors = 13
samp_neighbors = 11
mode = 'average'

# I think we already have it.
sc.pp.neighbors(thisq, n_neighbors=num_neighbors)

scores = []
for celli in range(0, len(thisq)):
    scores.append(sc_score_list(thisq, celli,
                                  noise_trials, num_neighbors, samp_neighbors, set_list,
                                  mode, False)


                  )

print(scores)
#thisq.obs['module_A'] = scores

#sc.pl.violin(thisq, keys='module_A', groupby='diagnosis')