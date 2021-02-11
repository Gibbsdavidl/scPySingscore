
import scanpy as sc
from scsingscore.scsingscore import score
import pandas as pd


def add_noise(df, n, noise_low, noise_high):
    df2 = df.copy()
    for i in range(0,n):
        df2['trial'+str(i)] = df * (1 + np.random.uniform(noise_low,noise_high,(df.shape)))
    return(df2)


def get_genecounts(adata, celli, count_threshold):
    avar = adata.var.copy()
    avar = avar.drop(['gene_ids', 'feature_types', 'genome', 'is_mito', 'is_ribo'], axis=1)
    gene_counts = adata.X[celli,:]
    if "toarray" in dir(gene_counts):
        avar['genecounts'] = gene_counts.toarray()[0]
    else:
        avar['genecounts'] = gene_counts.tolist()[0]
    avar_cut = avar[avar['genecounts'] > count_threshold]  ### MANY duplicate gene counts (i.e. 1)
    return(avar_cut)


def get_genecounts_csamp(adata, csamp, count_threshold):
    # get gene counts for a neighborhood sample
    avar = adata.var.copy()
    avar = avar.drop(['gene_ids', 'feature_types', 'genome', 'is_mito', 'is_ribo'], axis=1)
    gene_counts = adata.X[csamp,:]
    #if "toarray" in dir(gene_counts):
    #    avar['genecounts'] = gene_counts.toarray()[0]
    #else:
    #    avar['genecounts'] = gene_counts.tolist()[0]
    #avar_cut = avar[avar['genecounts'] > count_threshold]  ### MANY duplicate gene counts (i.e. 1)
    return(avar_cut)



def get_conn_dist(q, celli, nn):
    # q is an adata
    # celli is cell i in q
    # nn is the number of neighbors
    jdx = np.where(q.obsp['distances'][celli].todense() > 0)[1]
    y = [q.obsp['distances'][celli].todense().tolist()[0][ i ] for i in jdx]
    x = [q.obsp['connectivities'][celli].todense().tolist()[0][ i ] for i in jdx]
    df = pd.DataFrame({'idx':jdx, 'conn':x, 'dist':y})
    sumconn = sum(df['conn'])
    df['prob'] = [x/sumconn for x in df['conn']]
    return(df)


# In[50]:


def getIndices(qx, gs):
    idx = []
    for gi in gs:
        res0 = np.where(qx.var.index == gi)[0]
        if len(res0) > 0:
            idx.append(res0[0])
    return(idx)

def simValueHigh(qx, gs):
    # first get the index to X for these genes
    idx = getIndices(qx, gs)
    qy = qx[:,idx[20]].X
    return(qy)


# In[51]:


def neighborhood_scores(adata, noise_trials, num_neighbors, samp_neighbors, gene_set, mode):

    # mode 'average' averaged the noise trials
    # mode 'nonoise' returns the non-noised score 
    
    res0 = []
    numcells = adata.obs.shape[0]
    sc.pp.neighbors(adata, n_neighbors = num_neighbors)

    for celli in range(0, numcells):

        # first we get the neighborhood cells
        cdf = get_conn_dist(adata, celli, num_neighbors)

        # then we sample a set of cells from them
        csamp = np.random.choice(a=cdf['idx'], size=samp_neighbors, replace=False, p=cdf['prob'])

        # gene counts for scoring needs to have genes on rows.
        gene_counts = adata.X[csamp,:]
        gene_mat = gene_counts.todense()
        gene_mat = gene_mat.transpose().tolist()

        # get index of genes with some counts across cells
        gdx = [i for i,j in enumerate(gene_mat) if sum(j) > 0.0]  ## could filter here

        # sum gene counts across cells
        gene_mat_sum = [sum(x) for x in gene_mat]
        df = pd.DataFrame(gene_mat_sum, index=adata.var.index)
        df = df.iloc[gdx,:]

        if mode == 'average':
            # add some noise to gene counts.. create a n numbers of examples
            df_noise = add_noise(df, noise_trials, 0.01, 0.49)
            # score the neighborhoods
            si = score(up_gene=gene_set, sample=df_noise, norm_method='standard', full_data=False)  # standard workin gbetter here than theoretical
            res0.append(sum(si['total_score'][1:])/ float(noise_trials))
        else:
            si = score
            res0.append(float(si['total_score'][0]))
            
    return(res0)


# In[52]:


def one_score(adata, celli, noise_trials, num_neighbors, samp_neighbors, gene_set, 
              mode='average', compute_neighbors=True):

    # mode 'average' averaged the noise trials
    # mode 'nonoise' returns the non-noised score 
    
    if num_neighbors == 0 and samp_neighbors > 0:
        print('fix parameters')
        return(False)
    
    numcells = adata.obs.shape[0]
    
    if compute_neighbors == True and num_neighbors > 0:
        sc.pp.neighbors(adata, n_neighbors = num_neighbors)

    if num_neighbors > 0:
        # first we get the neighborhood cells
        cdf = get_conn_dist(adata, celli, num_neighbors)
        # then we sample a set of cells from them
        csamp = np.random.choice(a=cdf['idx'], size=samp_neighbors, replace=False, p=cdf['prob'])
        # gene counts for scoring needs to have genes on rows.
        gene_counts = adata.X[csamp,:]
    else:
        # one cell
        gene_counts = adata.X[celli,:]
        
    gene_mat = gene_counts.todense()    
    gene_mat = gene_mat.transpose().tolist()
    
    # get index of genes with some counts across cells
    gdx = [i for i,j in enumerate(gene_mat) if sum(j) > 0.0]  ## could filter here
    # sum gene counts across cells
    gene_mat_sum = [sum(x) for x in gene_mat]
    df = pd.DataFrame(gene_mat_sum, index=adata.var.index)
    df = df.iloc[gdx,:]
        
    if mode == 'average' and noise_trials > 0:
        # add some noise to gene counts.. create a n numbers of examples
        df_noise = add_noise(df, noise_trials, 0.01, 0.99)
        # score the neighborhoods
        si = score(up_gene=gene_set, sample=df_noise, norm_method='standard', full_data=False)  # standard workin gbetter here than theoretical
    else:
        si = score(up_gene=gene_set, sample=df, norm_method='standard', full_data=False) 
        
    return(si)


