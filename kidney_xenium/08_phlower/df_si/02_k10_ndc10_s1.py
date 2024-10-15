#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from collections import Counter
from datetime import datetime

import phlower
import importlib
importlib.reload(phlower)
phlower.__version__


adata = pickle.load(open("save/phlower_k10_ndc10_s1.pickle", 'rb'))
adata.obs['group'] = adata.obs['annotation']
phlower.tl.construct_delaunay(adata, trunc_quantile = 0.9, trunc_times = 3, start_n=10, end_n=10, separate_ends_triangle=True, circle_quant=0.7, calc_layout=False)

dt1 = datetime.now()
print('L1 Datetime Start:', dt1)

phlower.tl.L1Norm_decomp(adata, L1_mode='sym', check_symmetric=False, isnorm = True)
dt2 = datetime.now()
print('L1 Datetime End:', dt2)
delta = dt2 - dt1
print('L1 Difference is:', delta)

phlower.tl.knee_eigen(adata)
phlower.tl.random_climb_knn(adata, n=10000, roots_ratio=0.0001, knn_edges_k=20)
phlower.tl.trajs_matrix(adata)

phlower.tl.trajs_clustering(adata, eps=0.5, oname_basis='')
Counter(adata.uns['trajs_clusters'])

phlower.tl.select_trajectory_clusters(adata, skewness_threshold=2, manual_rm_clusters=[1])

phlower.tl.harmonic_stream_tree(adata, trajs_clusters='trajs_clusters',  eigen_n=-1, trajs_use=10000, kde_sample_n=10000, pca_name='X_pca', trim_end=False, min_bin_number=10, cut_threshold=1, verbose=True)


pickle.dump(adata, open("save/phlower_k10_ndc10_s1_step2.pickle", 'wb'))

