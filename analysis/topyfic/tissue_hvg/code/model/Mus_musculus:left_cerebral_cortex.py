import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Mus_musculus:left_cerebral_cortex')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Mus_musculus:left_cerebral_cortex')

ENCSR720TUP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR720TUP.p")
ENCSR632EIK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR632EIK.p")
ENCSR133XJX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR133XJX.p")
ENCSR624WDN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR624WDN.p")
ENCSR467GGT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR467GGT.p")
ENCSR224RWH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR224RWH.p")
ENCSR421YXX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR421YXX.p")
ENCSR132MEW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR132MEW.p")
ENCSR490RCS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR490RCS.p")
ENCSR037LFG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR037LFG.p")
ENCSR827AOS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR827AOS.p")
ENCSR527LIA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR527LIA.p")
ENCSR084RQP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR084RQP.p")
ENCSR118FCJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR118FCJ.p")
ENCSR611QJH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR611QJH.p")
ENCSR248MFG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR248MFG.p")
ENCSR450IJU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR450IJU.p")
ENCSR254SLA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR254SLA.p")
ENCSR839BHY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR839BHY.p")
ENCSR793HMA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR793HMA.p")
ENCSR308DHB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR308DHB.p")
ENCSR883XZK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR883XZK.p")
ENCSR632SJM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR632SJM.p")
ENCSR220GBV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR220GBV.p")
ENCSR257BOU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR257BOU.p")
ENCSR386PPF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR386PPF.p")
ENCSR081SXM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR081SXM.p")
ENCSR659DFF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR659DFF.p")
ENCSR964OQW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR964OQW.p")
ENCSR544MGC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR544MGC.p")
ENCSR302VFA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR302VFA.p")
ENCSR132KMR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR132KMR.p")
ENCSR178PKG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR178PKG.p")
ENCSR250SCW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR250SCW.p")
ENCSR092RYV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR092RYV.p")
ENCSR134FSZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR134FSZ.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_hvg_norm/Mus_musculus:left_cerebral_cortex.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR720TUP,ENCSR632EIK,ENCSR133XJX,ENCSR624WDN,ENCSR467GGT,ENCSR224RWH,ENCSR421YXX,ENCSR132MEW,ENCSR490RCS,ENCSR037LFG,ENCSR827AOS,ENCSR527LIA,ENCSR084RQP,ENCSR118FCJ,ENCSR611QJH,ENCSR248MFG,ENCSR450IJU,ENCSR254SLA,ENCSR839BHY,ENCSR793HMA,ENCSR308DHB,ENCSR883XZK,ENCSR632SJM,ENCSR220GBV,ENCSR257BOU,ENCSR386PPF,ENCSR081SXM,ENCSR659DFF,ENCSR964OQW,ENCSR544MGC,ENCSR302VFA,ENCSR132KMR,ENCSR178PKG,ENCSR250SCW,ENCSR092RYV,ENCSR134FSZ], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

