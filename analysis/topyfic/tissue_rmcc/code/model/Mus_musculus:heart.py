import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Mus_musculus:heart')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Mus_musculus:heart')

ENCSR333SLJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR333SLJ.p")
ENCSR415HEP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR415HEP.p")
ENCSR568GGW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR568GGW.p")
ENCSR533IGQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR533IGQ.p")
ENCSR665HPR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR665HPR.p")
ENCSR417FCM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR417FCM.p")
ENCSR462HCD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR462HCD.p")
ENCSR644VVE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR644VVE.p")
ENCSR729BRE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR729BRE.p")
ENCSR266RTS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR266RTS.p")
ENCSR442NQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR442NQZ.p")
ENCSR440PMT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR440PMT.p")
ENCSR938YKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR938YKJ.p")
ENCSR041KFB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR041KFB.p")
ENCSR906KCJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR906KCJ.p")
ENCSR764CYB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR764CYB.p")
ENCSR594LDC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR594LDC.p")
ENCSR860LQG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR860LQG.p")
ENCSR303EIR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR303EIR.p")
ENCSR617CNE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR617CNE.p")
ENCSR305SPL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR305SPL.p")
ENCSR147TEP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR147TEP.p")
ENCSR122DEZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR122DEZ.p")
ENCSR247PQY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR247PQY.p")
ENCSR868XSG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR868XSG.p")
ENCSR975YAQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR975YAQ.p")
ENCSR125HRM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR125HRM.p")
ENCSR277UKO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR277UKO.p")
ENCSR237OCP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR237OCP.p")
ENCSR688APT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR688APT.p")
ENCSR559AAG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR559AAG.p")
ENCSR337DKH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR337DKH.p")
ENCSR809WJL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR809WJL.p")
ENCSR573UKV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR573UKV.p")
ENCSR743OKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR743OKJ.p")
ENCSR455YDL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR455YDL.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_rmcc_norm/Mus_musculus:heart.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR333SLJ,ENCSR415HEP,ENCSR568GGW,ENCSR533IGQ,ENCSR665HPR,ENCSR417FCM,ENCSR462HCD,ENCSR644VVE,ENCSR729BRE,ENCSR266RTS,ENCSR442NQZ,ENCSR440PMT,ENCSR938YKJ,ENCSR041KFB,ENCSR906KCJ,ENCSR764CYB,ENCSR594LDC,ENCSR860LQG,ENCSR303EIR,ENCSR617CNE,ENCSR305SPL,ENCSR147TEP,ENCSR122DEZ,ENCSR247PQY,ENCSR868XSG,ENCSR975YAQ,ENCSR125HRM,ENCSR277UKO,ENCSR237OCP,ENCSR688APT,ENCSR559AAG,ENCSR337DKH,ENCSR809WJL,ENCSR573UKV,ENCSR743OKJ,ENCSR455YDL], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

