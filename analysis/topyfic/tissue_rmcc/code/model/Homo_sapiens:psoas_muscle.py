import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Homo_sapiens:psoas_muscle')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Homo_sapiens:psoas_muscle')

ENCSR471YTA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR471YTA.p")
ENCSR203THX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR203THX.p")
ENCSR005TCL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR005TCL.p")
ENCSR281RIH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR281RIH.p")
ENCSR641MWO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR641MWO.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_rmcc_norm/Homo_sapiens:psoas_muscle.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR471YTA,ENCSR203THX,ENCSR005TCL,ENCSR281RIH,ENCSR641MWO], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

