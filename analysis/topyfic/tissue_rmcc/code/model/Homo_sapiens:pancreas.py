import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Homo_sapiens:pancreas')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Homo_sapiens:pancreas')

ENCSR472RRP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR472RRP.p")
ENCSR478FQR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR478FQR.p")
ENCSR261BXB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR261BXB.p")
ENCSR281NBH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR281NBH.p")
ENCSR684KYI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR684KYI.p")
ENCSR445LCA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR445LCA.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_rmcc_norm/Homo_sapiens:pancreas.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR472RRP,ENCSR478FQR,ENCSR261BXB,ENCSR281NBH,ENCSR684KYI,ENCSR445LCA], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

