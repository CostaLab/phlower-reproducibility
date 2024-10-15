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

import pickle
pickle.dump(adata, open(f"save/phlower_d12_18.pickle", 'wb'))
