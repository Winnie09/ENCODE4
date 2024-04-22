import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Homo_sapiens:left_lung')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Homo_sapiens:left_lung')

ENCSR966DDY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR966DDY.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_hvg_norm/Homo_sapiens:left_lung.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR966DDY], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

