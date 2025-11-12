import os
import numpy as np
from utils.plot import plot
from scDML import scDMLModel

def scdml(adata, out_path, batch, celltype):
    out_path = out_path + "scdml"

    adata.obs["batch_orig"] = adata.obs[batch].copy()

    model = scDMLModel()
    adata = model.preprocess(adata, batch_key = batch)
    model.integrate(adata, 
                    batch_key = batch, 
                    ncluster_list = [15], 
                    merge_rule="rule2", 
                    expect_num_cluster = 15, 
                    mode = "unsupervised")
    
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.save(f"{out_path}_emb.npy", adata.obsm['X_emb'])
    
    plot(adata, "X_emb", 'batch_orig' , celltype, out_path)
    return adata