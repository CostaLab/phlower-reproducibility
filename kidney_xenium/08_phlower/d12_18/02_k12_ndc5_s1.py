#!/usr/bin/env python
# coding: utf-8
from datetime import datetime
import phlower
import importlib
importlib.reload(phlower)
phlower.__version__


adata = pickle.load(open("save/phlower_k12_ndc5_s1.pickle", 'rb'))
adata.obs['group'] = adata.obs['annotation']
phlower.tl.construct_delaunay(adata, trunc_quantile = 0.9, trunc_times = 3, start_n=10, end_n=10, separate_ends_triangle=True, circle_quant=0.7, calc_layout=False)

fig, ax = plt.subplots(1,2, figsize=(20,8))
phlower.pl.nxdraw_group(adata, node_size=20, show_edges=True, ax=ax[0], show_legend=False, labelstyle='box', labelsize=8)
phlower.pl.nxdraw_group(adata, graph_name='X_pca_ddhodge_g_triangulation_circle', node_size=20, show_edges=True, ax=ax[1], show_legend=True, labelstyle='box', labelsize=8)



dt1 = datetime.now()
print('L1 Datetime Start:', dt1)

phlower.tl.L1Norm_decomp(adata, L1_mode='sym', check_symmetric=False, isnorm = True)
dt2 = datetime.now()
print('L1 Datetime End:', dt2)
delta = dt2 - dt1
print('L1 Difference is:', delta)

phlower.tl.knee_eigen(adata)
phlower.tl.random_climb_knn(adata, n=10000, roots_ratio=0.0001, knn_edges_k=15)
phlower.tl.trajs_matrix(adata)
phlower.tl.trajs_clustering(adata, eps=0.5, oname_basis='')
Counter(adata.uns['trajs_clusters'])
phlower.tl.select_trajectory_clusters(adata, skewness_threshold=2)
phlower.tl.harmonic_stream_tree(adata, trajs_clusters='trajs_clusters',  eigen_n=-1, trajs_use=10000, kde_sample_n=10000, pca_name='X_pca', trim_end=False, min_bin_number=20, cut_threshold=1, verbose=True)
pickle.dump(adata, open("save/phlower_k12_ndc5_s1_step2.pickle", 'wb'))


