import os
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
import logging
import anndata2ri
import rpy2
import rpy2.rinterface_lib.callbacks
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)  # Ignore R warning messages
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
import warnings
warnings.filterwarnings('ignore')

ro.r.source('utils/batchKL.R')
ro.r.source('utils/calLISI.R')

from time import time
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import pair_confusion_matrix

# sklearn ari bug
def ari(labels_true, labels_pred): 
    (tn, fp), (fn, tp) = pair_confusion_matrix(labels_true, labels_pred)
    tn=int(tn)
    tp=int(tp)
    fp=int(fp)
    fn=int(fn)

    if fn == 0 and fp == 0:
        return 1.0

    return 2. * (tp * tn - fn * fp) / ((tp + fn) * (fn + tn) + (tp + fp) * (fp + tn))


### kBET
def kBET_single(matrix, batch, type_ = None, k0 = 20, knn=None, subsample=0.5, heuristic=True, verbose=False):
    """
    params:
        matrix: expression matrix (currently a PCA matrix)
        batch: series or list of batch assignments
        subsample: fraction to be subsampled. No subsampling if `subsample=None`
    returns:
        kBET p-value
    """
    ro.r("library(kBET)")

    if verbose:
        print("importing expression matrix")
    ro.globalenv['data_mtrx'] = matrix
    ro.globalenv['batch'] = batch
    ro.globalenv['k0'] = k0

    try:
        batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, k0=k0, plot=FALSE, do.pca=FALSE)")
        result = ro.r("mean(batch.estimate$stats$kBET.observed, na.rm=T)")
        return result
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        return np.nan


### BatchKL
def BatchKL(adata_integrated):
    embedding = adata_integrated.obsm["X_emb"]

    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        meta_data = ro.conversion.py2rpy(adata_integrated.obs)
        KL = ro.r.BatchKL(meta_data, embedding, n_cells=100, batch="BATCH")
    
    print("BatchKL =", KL)
    return KL


### LISI
def LISI(adata_integrated):
    embedding = adata_integrated.obsm["X_emb"]

    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        meta_data = ro.conversion.py2rpy(adata_integrated.obs)
        lisi = ro.r.CalLISI(embedding, meta_data)
    
    print("cLISI =", lisi[0])
    print("iLISI =", lisi[1])
    return lisi


### Silhouette score
def silhouette(adata, group_key, embed, metric='euclidean', scale=True):
    if embed not in adata.obsm.keys():
        raise KeyError(f'{embed} not in adata.obsm')
    asw = sklearn.metrics.silhouette_score(
        X=adata.obsm[embed],
        labels=adata.obs[group_key],
        metric=metric
    )
    if scale:
        asw = (asw + 1)/2
    return asw


def silhouette_batch(adata, batch_key, group_key, embed, metric='euclidean', verbose=True, scale=True):
    if embed not in adata.obsm.keys():
        raise KeyError(f'{embed} not in adata.obsm')
    
    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])
    
    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()
        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue
        sil_per_group = sklearn.metrics.silhouette_samples(
            adata_group.obsm[embed], adata_group.obs[batch_key], metric=metric
        )
        sil_per_group = [abs(i) for i in sil_per_group]
        if scale:
            sil_per_group = [1 - i for i in sil_per_group]
        d = pd.DataFrame({'group': [group]*len(sil_per_group), 'silhouette_score': sil_per_group})
        sil_all = pd.concat([sil_all, d])
    
    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    
    if verbose:
        print(f'mean silhouette per cell: {sil_means}')
    return sil_all, sil_means


# Find optimal resolution
def find_resolution(adata_, n_clusters, random):
    adata = adata_.copy()
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1000.]
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        sc.tl.louvain(adata, resolution=current_res, random_state=random)
        labels = adata.obs['louvain']
        obtained_clusters = len(np.unique(labels))

        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        iteration += 1
    return current_res


def calulate_ari_nmi(adata_integrated, n_cluster=4):
    sc.pp.neighbors(adata_integrated, random_state=0)
    reso = find_resolution(adata_integrated, n_cluster, 0)
    sc.tl.louvain(adata_integrated, reso, random_state=0)
    sc.tl.umap(adata_integrated)

    if adata_integrated.X.shape[1] == 2:
        adata_integrated.obsm["X_emb"] = adata_integrated.X

    ARI = ari(adata_integrated.obs["celltype"].astype(str), adata_integrated.obs["louvain"])
    NMI = normalized_mutual_info_score(
        adata_integrated.obs["celltype"].astype(str), adata_integrated.obs["louvain"]
    )

    print(f"louvain clustering result(resolution={reso}): n_cluster={n_cluster}")
    print("ARI:", ARI)
    print("NMI:", NMI)
    return ARI, NMI


def evaluate_dataset(adata_integrated, method="louvain1.0", n_celltype=10):
    print(f"method={method}, n_celltype={n_celltype}")
    ARI, NMI = calulate_ari_nmi(adata_integrated, n_cluster=n_celltype)

    KL = BatchKL(adata_integrated)
    lisi = LISI(adata_integrated)

    sil_global = silhouette(adata_integrated, group_key="celltype", embed="X_emb", metric='euclidean')
    _, sil_clus = silhouette_batch(
        adata_integrated, batch_key="BATCH", group_key="celltype", embed="X_emb",
        metric='euclidean', verbose=False
    )
    sil_clus = sil_clus['silhouette_score'].mean()

    results = {
        'ARI': np.round(ARI, 3),
        'NMI': np.round(NMI, 3),
        'ASW_label': np.round(sil_global, 3),
        'ASW_label/batch': np.round(sil_clus, 3),
        'BatchKL': np.round(KL[0], 3),
        'cLISI': np.round(lisi[0], 3),
        'iLISI': np.round(lisi[1], 3)
    }

    result = pd.DataFrame.from_dict(results, orient='index')
    result.columns = [method]
    print("All metrics calculation done.")
    return adata_integrated, result
