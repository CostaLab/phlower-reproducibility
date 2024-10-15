#!/usr/bin/env python
# coding: utf-8

# # Kidney

import os
import phlower
import numpy as np
import pandas as pd
import networkx as nx
from itertools import chain
import scipy
import pickle
import sklearn
import pydot
import scanpy as sc
import seaborn as sns
import colorcet as cc
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from datetime import datetime
from collections import Counter
import networkx as nx
import phlower





knn_k = 6
ndc_n = 40
s_n = 4



import logging

if not os.path.exists("logs"):
    os.mkdir("logs")
logger = logging.getLogger()
fhandler = logging.FileHandler(filename=f'logs/i_{knn_k}_{datetime.now()}.log'.replace(" ", "_"), mode='a')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

logger.info(f'begining: {datetime.now()}')


# ### load kidney anndata with MOJITOO reduction and clustering.


adata = phlower.dataset.kidney()


dt1 = datetime.now()
print('ddhodge Datetime Start:', dt1)
logger.info(f'ddhodge Datetime Start: {dt1}')

phlower.ext.ddhodge(adata, basis='X_pca', npc=30, roots=adata.obs['group']==6, k=knn_k, ndc=ndc_n,s=s_n)

dt2 = datetime.now()

print('ddhodge Datetime End:', dt2)
logger.info(f'ddhodge Datetime End: {dt2}')
delta = dt2 - dt1

print('ddhodge Difference is:', delta)
logger.info(f'ddhodge Difference is: {delta}')


# ### show the pseudotime on the neato layout




phlower.pl.nxdraw_score(adata, color='u',node_size=10)


# ### show the clusters




_, ax = plt.subplots(1,1, )
phlower.pl.nxdraw_group(adata,node_size=10, show_edges=False, ax=ax, legend_col=2)
ax.legend(ncol=2, markerscale=1, loc="center left", bbox_to_anchor=(1, 0.5))


# ### show the days




dic =  {"S1":"Day7_18", "S2":"Day7_12", "S3":"Day7_5", "S4":"Day7_0"}
adata.obs['time'] = [dic[i] for i in adata.obs['sample']]
phlower.pl.nxdraw_group(adata,node_size=10, group_name='time',show_edges=False, show_legend=True, legend_col=1, markerscale=3, label=False)


# ### connect start state cells and terminal state cells




phlower.tl.construct_delaunay(adata, trunc_quantile = 0.9, trunc_times = 3, start_n=10, end_n=10, separate_ends_triangle=True, circle_quant=0.1, calc_layout=False)


# ### Perform hodge laplacian decomposition




from datetime import datetime

dt1 = datetime.now()
print('L1 Datetime Start:', dt1)
logging.info(f'L1 Datetime Start: {dt1}')

phlower.tl.L1Norm_decomp(adata, L1_mode='sym', isnorm = True)
dt2 = datetime.now()
print('L1 Datetime End:', dt2)
logger.info(f'L1 Datetime End: {dt2}')
delta = dt2 - dt1
print('L1 Difference is:', delta)
logger.info(f'L1 Difference: {delta}')


# ### check the number of harmonics




_, ax = plt.subplots(1,2, figsize=(9,3))
phlower.pl.plot_eigen_line(adata, n_eig=10,ax=ax[0])
phlower.pl.plot_eigen_line(adata, n_eig=100, step_size=10, ax=ax[1])





phlower.tl.knee_eigen(adata)


# ### preference random walk to sample trajectories




phlower.tl.random_climb_knn(adata, n=10000)


# ### display 2 sampled trjaectoris




from scipy import constants
fig_width = 9
fig, axs = plt.subplots(1, 4,figsize=(fig_width*1.5, fig_width*0.5 / constants.golden))
phlower.pl.plot_traj(adata, trajectory=adata.uns['knn_trajs'][0], colorid=0, node_size=1, ax=axs[0])
phlower.pl.plot_traj(adata, trajectory=adata.uns['knn_trajs'][0], colorid=0, node_size=0, ax=axs[1])
phlower.pl.plot_traj(adata,  trajectory=adata.uns['knn_trajs'][1], colorid=0, node_size=1, ax=axs[2])
phlower.pl.plot_traj(adata,  trajectory=adata.uns['knn_trajs'][1], colorid=0, node_size=0, ax=axs[3])


# ### project trajectories to harmonic space




phlower.tl.trajs_matrix(adata)


# ### clustering trajectory paths in t-map space




phlower.tl.trajs_clustering(adata, eps=0.5, oname_basis='')


# ### clean bad clusters




phlower.tl.select_trajectory_clusters(adata, manual_rm_clusters=[2], verbose=False)
phlower.tl.unique_trajectory_clusters(adata, group_name='group', verbose=False)


# ### show trajectory clusters in umap (derived from t-map space)




_, ax = plt.subplots(1,1, figsize=(7, 4))
phlower.pl.plot_trajs_embedding(adata, clusters="trajs_clusters", ax=ax)


# ### show sampled trajectory paths in neato layout space




fig, ax = phlower.pl.plot_density_grid(adata,  figsize=(10, 7), cluster_name='trajs_clusters', bg_alpha=0.1, return_fig=True, sample_n=1000)





fig, ax = plt.subplots(1,2, figsize=(10,4))
phlower.pl.nxdraw_group(adata, node_size=20, show_edges=True, ax=ax[0], show_legend=False, labelstyle='text', labelsize=8)
phlower.pl.nxdraw_group(adata, node_size=20, show_edges=True, ax=ax[1], show_legend=True, labelstyle='box', labelsize=8, legend_col=2)


# ### annotate the branches




anno_dic = {
    0:  "Neuron-2",
    1:  "Tubular",
    3:  "Stromal-3",
    5:  "Stromal-4",
    6:  "Podocytes",
    7:  "Stromal-2",
    9:  "Muscle",
    11: "Stromal-1",
    10: "Neuron-1",
}
adata.uns['annotation']= [anno_dic.get(int(i), i) for i in adata.uns['trajs_clusters']]


# ### create stream tree




phlower.tl.harmonic_stream_tree(adata,
                                trajs_clusters='annotation',
                                retain_clusters=list(set(anno_dic.values())),
                                min_bin_number=20,
                                cut_threshold=1.5,
                                verbose=True)





adata.uns['workdir'] = ''
adata.obs['group_str'] = [str(i) for i in adata.obs['group']]


# ### plot stream tree




phlower.ext.plot_stream_sc(adata, fig_size=(8,5), color=['group_str'],show_legend=False, dist_scale=1,s=10)


# ### show trajectory paths in umap and in ct-map space




figs = []
fig, ax = plt.subplots(1,2, figsize=(6,3))
phlower.pl.plot_trajs_embedding(adata, clusters="annotation",show_legend=False, labelsize=5, node_size=10, ax=ax[0])
phlower.pl.plot_trajectory_harmonic_lines(adata, clusters="annotation", sample_ratio=0.05, show_legend=True, ax=ax[1])


# ### show trajctory paths in t-map space




phlower.pl.plot_trajectory_harmonic_points_3d(adata, clusters="annotation", sample_ratio=.05, show_legend=True, fig_path=None)


# ### show trajctory paths in ct-map space




phlower.pl.plot_trajectory_harmonic_lines_3d(adata, clusters="annotation", sample_ratio=.05, show_legend=True, fig_path=None)


# ### save adata with networkx object using pickle




import pickle
import os
from datetime import date

if not os.path.exists("save"):
    os.mkdir('save')

today_ = date.today().strftime("%Y-%m-%d")
pickle.dump(adata, open(f"save/kidney_{today_}.pickle", 'wb'))


# ### show session info




import session_info
session_info.show()

