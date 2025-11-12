import os
import scvi
import numpy as np
from utils.plot import plot

def scVi(adata, out_path, batch, celltype):
    out_path = out_path + "scvi"

    # work on a copy and return it so main can evaluate metrics
    adata = adata.copy()

    scvi.model.SCVI.setup_anndata(adata, batch_key = batch)
    model = scvi.model.SCVI(adata)
    model.train()
    adata.obsm["X_scVI"] = model.get_latent_representation()

    adata.obsm['X_emb'] = adata.obsm['X_scVI']

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.save(f"{out_path}_emb.npy", adata.obsm['X_emb'])

    plot(adata, "X_scVI", batch , celltype, out_path)
    return adata
