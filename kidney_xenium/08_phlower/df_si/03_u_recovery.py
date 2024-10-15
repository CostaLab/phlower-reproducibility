#!/usr/bin/env python
# coding: utf-8

# ### recover cells not included in phlower

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from collections import Counter, defaultdict
from datetime import datetime
import networkx as nx
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm


import phlower
import importlib
importlib.reload(phlower)
phlower.__version__



print(datetime.now(), "loading adata...")
adata = pickle.load(open("save/phlower_k10_ndc10_s1_step2.pickle", 'rb'))
print(datetime.now(), "loaded.")

print(datetime.now(), "loading bdata...")
bdata = pickle.load(open('save/phlower_df_si.pickle', 'rb'))
print(datetime.now(), "loaded.")



def branch_nodes(bfrom='4_18', tree='fate_tree', graph_name = None):
    graph_name = adata.uns["graph_basis"] + "_triangulation_circle"
    elist = np.array([(x[0], x[1]) for x in adata.uns[graph_name].edges()])
    edict = {i:v for i, v in enumerate(elist)}
    nodes = list(nx.dfs_tree(adata.uns[tree], bfrom))
    print(nodes)
    node_dict = defaultdict(int)
    for i in nodes:
        ecount = adata.uns[tree].nodes[i]['ecount']
        for k,v in ecount:
            n1,n2 = edict[k]
            node_dict[n1] += v
            node_dict[n2] += v

    return node_dict



def a2b_data_index(adata, bdata, indices):
    """
    adata is a subset of bdata
    given adata indices, return bdata indices
    """
    bc = adata.obs_names.take(indices)
    bidxs = [list(bdata.obs_names).index(x) for x in bc]
    return {aidx:bidx for aidx,bidx in zip(indices, bidxs)}



def celltype(adata, indices, name='annotation'):
    """
    return celltype of adata indices
    """
    return adata.obs[name].cat.categories[adata.obs[name].cat.codes.take(indices)]

def idxstocelltypes(adata, indices):
    return celltype(adata, indices)

## for each cell type, find its neighbors in the subset, and get the average pseudotime.
knn=10
dm='X_pca'
br = 'u'
decay=np.sqrt
#clip=10
knn_k = 3000
name = 'annotation'

node_dict = {i:e for i, e in enumerate(adata.obs['u'])}


node_dict = {k:int(decay(v)) for k,v in node_dict.items()}
pseudo_dict = {i: adata.obs[br].iloc[i] for i in node_dict}
print(datetime.now(), "knn...")
nbrs = NearestNeighbors(n_neighbors=3000, algorithm='ball_tree', n_jobs=100).fit(bdata.obsm[dm])
print(datetime.now(), "knn done.")

##bdata.obs['annotation'] to category
bdata.obs['annotation'] = bdata.obs['annotation'].astype('category')

if True:
    ddd = a2b_data_index(adata,bdata, list(node_dict.keys()))
    rddd = {v:k for k,v in ddd.items()}
    indist, inind = nbrs.kneighbors(bdata.obsm[dm])
    pickle.dump({"indist":indist, "inind":inind, 'a2b_index':ddd, "r_a2b_index":rddd}, open("save/indist_inind_df_si.pickle", 'wb'))
else:
    d = pickle.load(open("save/indist_inind_df_si.pickle", 'rb'))
    indist, inind, ddd, rddd = d['indist'], d['inind'], d['a2b_index'], d['r_a2b_index']

rest = set(range(bdata.n_obs)) - (set(ddd.values()))
is_neigbor_set = set()
pseudotime_dict = {} ##
print(datetime.now(), "psuedo calc...")
for i in tqdm(range(bdata.n_obs)):
    if i in rddd.keys(): ## already in subset
        pseudotime_dict[i] = pseudo_dict[rddd[i]]
        is_neigbor_set.add(i)
        continue
    ct = celltype(bdata, [i])[0]
    nns = inind[i, 1:(knn_k+1)]
    nns_ct = idxstocelltypes(bdata, nns)
    nns =nns[np.where(nns_ct == ct)[0]]
    if len(nns) == 0:
        print("no ct neighbors:",  nns)
        break
    #print(len(nns))

    u = set(nns) & set(rddd.keys())
    u = list(u)[:30]
    length = len(u)
    if length > 0:
        aidxs = [rddd[i] for i in u]
        is_neigbor_set.add(i)
        #print(aidxs)
        pseudotime_dict[i] = np.mean([pseudo_dict[i] for i in aidxs])
    else:
        print("no neighbors:", i)

print(datetime.now(), "psuedo done.")
bdata.obs[f'{name}_pseudo'] = -1
bdata.obs[f'{name}_pseudo'][list(pseudotime_dict.keys())] = list(pseudotime_dict.values())


import pickle
pickle.dump(bdata, open("save/03_bdata_u_recovery.pickle", "wb"))

bdata.obs.to_csv("save/03_bdata_u_recovery.csv")
