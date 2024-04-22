import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Mus_musculus:layer_of_hippocampus')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Mus_musculus:layer_of_hippocampus')

ENCSR944XPC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR944XPC.p")
ENCSR995ZKH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR995ZKH.p")
ENCSR539XHZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR539XHZ.p")
ENCSR987SOK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR987SOK.p")
ENCSR835JFN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR835JFN.p")
ENCSR968RPR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR968RPR.p")
ENCSR742NWO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR742NWO.p")
ENCSR075HVN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR075HVN.p")
ENCSR079GGS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR079GGS.p")
ENCSR846QIZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR846QIZ.p")
ENCSR483SSM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR483SSM.p")
ENCSR642QUS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR642QUS.p")
ENCSR569LJC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR569LJC.p")
ENCSR399EUZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR399EUZ.p")
ENCSR967YET=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR967YET.p")
ENCSR940VCU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR940VCU.p")
ENCSR203XJO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR203XJO.p")
ENCSR435CMT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR435CMT.p")
ENCSR545URU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR545URU.p")
ENCSR257VQO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR257VQO.p")
ENCSR127LUZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR127LUZ.p")
ENCSR859XOI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR859XOI.p")
ENCSR736JKX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR736JKX.p")
ENCSR499ZYZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR499ZYZ.p")
ENCSR983TAK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR983TAK.p")
ENCSR613VKV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR613VKV.p")
ENCSR767FKE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR767FKE.p")
ENCSR824CFI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR824CFI.p")
ENCSR329JYG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR329JYG.p")
ENCSR079CKL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR079CKL.p")
ENCSR888TQC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR888TQC.p")
ENCSR979MRY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR979MRY.p")
ENCSR973WCI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR973WCI.p")
ENCSR247QPJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR247QPJ.p")
ENCSR730JFZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR730JFZ.p")
ENCSR167HUO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR167HUO.p")
ENCSR869HKZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR869HKZ.p")
ENCSR926QCE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR926QCE.p")
ENCSR692LZA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR692LZA.p")
ENCSR320EFC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR320EFC.p")
ENCSR232VMC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR232VMC.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_rmcc_norm/Mus_musculus:layer_of_hippocampus.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR944XPC,ENCSR995ZKH,ENCSR539XHZ,ENCSR987SOK,ENCSR835JFN,ENCSR968RPR,ENCSR742NWO,ENCSR075HVN,ENCSR079GGS,ENCSR846QIZ,ENCSR483SSM,ENCSR642QUS,ENCSR569LJC,ENCSR399EUZ,ENCSR967YET,ENCSR940VCU,ENCSR203XJO,ENCSR435CMT,ENCSR545URU,ENCSR257VQO,ENCSR127LUZ,ENCSR859XOI,ENCSR736JKX,ENCSR499ZYZ,ENCSR983TAK,ENCSR613VKV,ENCSR767FKE,ENCSR824CFI,ENCSR329JYG,ENCSR079CKL,ENCSR888TQC,ENCSR979MRY,ENCSR973WCI,ENCSR247QPJ,ENCSR730JFZ,ENCSR167HUO,ENCSR869HKZ,ENCSR926QCE,ENCSR692LZA,ENCSR320EFC,ENCSR232VMC], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

