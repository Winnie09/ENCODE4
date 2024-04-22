import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Mus_musculus:layer_of_hippocampus')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Mus_musculus:layer_of_hippocampus')

ENCSR944XPC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR944XPC.p")
ENCSR995ZKH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR995ZKH.p")
ENCSR539XHZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR539XHZ.p")
ENCSR987SOK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR987SOK.p")
ENCSR835JFN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR835JFN.p")
ENCSR968RPR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR968RPR.p")
ENCSR742NWO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR742NWO.p")
ENCSR075HVN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR075HVN.p")
ENCSR079GGS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR079GGS.p")
ENCSR846QIZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR846QIZ.p")
ENCSR483SSM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR483SSM.p")
ENCSR642QUS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR642QUS.p")
ENCSR569LJC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR569LJC.p")
ENCSR399EUZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR399EUZ.p")
ENCSR967YET=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR967YET.p")
ENCSR940VCU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR940VCU.p")
ENCSR203XJO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR203XJO.p")
ENCSR435CMT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR435CMT.p")
ENCSR545URU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR545URU.p")
ENCSR257VQO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR257VQO.p")
ENCSR127LUZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR127LUZ.p")
ENCSR859XOI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR859XOI.p")
ENCSR736JKX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR736JKX.p")
ENCSR499ZYZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR499ZYZ.p")
ENCSR983TAK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR983TAK.p")
ENCSR613VKV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR613VKV.p")
ENCSR767FKE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR767FKE.p")
ENCSR824CFI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR824CFI.p")
ENCSR329JYG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR329JYG.p")
ENCSR079CKL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR079CKL.p")
ENCSR888TQC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR888TQC.p")
ENCSR979MRY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR979MRY.p")
ENCSR973WCI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR973WCI.p")
ENCSR247QPJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR247QPJ.p")
ENCSR730JFZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR730JFZ.p")
ENCSR167HUO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR167HUO.p")
ENCSR869HKZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR869HKZ.p")
ENCSR926QCE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR926QCE.p")
ENCSR692LZA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR692LZA.p")
ENCSR320EFC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR320EFC.p")
ENCSR232VMC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR232VMC.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_hvg_norm/Mus_musculus:layer_of_hippocampus.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR944XPC,ENCSR995ZKH,ENCSR539XHZ,ENCSR987SOK,ENCSR835JFN,ENCSR968RPR,ENCSR742NWO,ENCSR075HVN,ENCSR079GGS,ENCSR846QIZ,ENCSR483SSM,ENCSR642QUS,ENCSR569LJC,ENCSR399EUZ,ENCSR967YET,ENCSR940VCU,ENCSR203XJO,ENCSR435CMT,ENCSR545URU,ENCSR257VQO,ENCSR127LUZ,ENCSR859XOI,ENCSR736JKX,ENCSR499ZYZ,ENCSR983TAK,ENCSR613VKV,ENCSR767FKE,ENCSR824CFI,ENCSR329JYG,ENCSR079CKL,ENCSR888TQC,ENCSR979MRY,ENCSR973WCI,ENCSR247QPJ,ENCSR730JFZ,ENCSR167HUO,ENCSR869HKZ,ENCSR926QCE,ENCSR692LZA,ENCSR320EFC,ENCSR232VMC], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

