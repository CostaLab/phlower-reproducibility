#!/usr/bin/env python
# coding: utf-8
import scipy.spatial
import sys
import getopt
import phlower
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

k=9
ndc=40
s=2

try:
    options,args = getopt.getopt(sys.argv[1:], "k:d:s:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for key, val in options:
    if key in "-k":
        k = val
        print("k: ", val)
    elif key in "-d":
        ndc = val
        print("ndc: ", val)
    elif key in "-s":
        s = val
        print("s: ", val)


print(f"k{k}_ndc{ndc}_s{s}")

meta = pd.read_csv("../../save/07_all_353_cleaned_subset_meta.csv", index_col=0)


meta = meta.loc[meta['sample'].isin(["iPS353_d12", "iPS353_d18"]), :]
meta.head(2)



from collections import Counter
Counter(meta['annotation'])



pcs = pd.read_csv("../../save/07_all_353_cleaned_subset_pca.csv", index_col=0) 
pcs = pcs.loc[meta.index, :]
pcs.head(2)




umap = pd.read_csv("../../save/07_all_353_cleaned_subset_umap.csv", index_col=0)
umap = umap.loc[meta.index,:]
umap.head(2)



rna_counts = sc.read_10x_h5('../../save/07_all_353_cleaned_pca_annotation_subset.h5')



adata = rna_counts[rna_counts.obs_names.isin(meta.index)]
adata

del rna_counts


adata.obsm['X_pca'] = np.array(pcs)
adata.obsm['X_umap'] = np.array(umap)
adata.obs['sample'] = meta['sample']
adata.obs['annotation'] = meta['annotation']
adata.obs['seurat_clusters'] = meta['seurat_clusters']

del pcs, umap, meta

sc.pl.umap(adata, color='annotation', legend_loc='on data')


def subset_celltype(adata, slot='annotation', minn=2000, ratio=0.2):
    np.random.seed(2024)
    if ratio >=1:
        return adata
    idx_list = []
    for ct in Counter(adata.obs[slot]).keys():
        if len(np.where(adata.obs[slot]==ct)[0]) > minn:
            aidx = np.where(adata.obs[slot] == ct)[0]
            idx = np.random.choice(aidx, max(int(len(aidx)*ratio), 1), replace=False)
        else:
            idx = np.where(adata.obs[slot] == ct)[0]
        idx_list.extend(idx)
    return adata[idx_list, :]
#endf subset_celltype


bdata = subset_celltype(adata, ratio=0.2)
bdata, Counter(bdata.obs['annotation'])

del adata

dt1 = datetime.now() 
print(dt1, flush=True)
phlower.ext.ddhodge(bdata, basis='X_pca',
                    npc=50, roots=bdata.obs['annotation']=='Mesoderm',
                    k=int(k), ndc=int(ndc),s=int(s),
                    lstsq_method='cholesky',
                    verbose=True)
dt2 = datetime.now() 
print(dt2, flush=True)


fig, ax = plt.subplots(2,1, figsize=[4, 8])
phlower.pl.nxdraw_group(bdata, group_name='annotation', layout_name='X_pca_ddhodge_g', node_size=5, show_edges=False, labelstyle='text', labelsize=6, ax=ax[0])
phlower.pl.nxdraw_score(bdata, layout_name='X_pca_ddhodge_g', node_size=5, ax=ax[1])
fig.savefig(f"viz/layout_k{k}_ndc{ndc}_s{s}.pdf", bbox_inches = 'tight')


import pickle
pickle.dump(bdata, open(f"save/phlower_k{k}_ndc{ndc}_s{s}.pickle", 'wb'))
