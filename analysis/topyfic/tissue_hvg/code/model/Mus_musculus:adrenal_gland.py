import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Mus_musculus:adrenal_gland')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_hvg/res/Mus_musculus:adrenal_gland')

ENCSR402PHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR402PHR.p")
ENCSR157YJO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR157YJO.p")
ENCSR021YYL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR021YYL.p")
ENCSR988AMY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR988AMY.p")
ENCSR908CQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR908CQZ.p")
ENCSR749GDE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR749GDE.p")
ENCSR224OUG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR224OUG.p")
ENCSR356VJZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR356VJZ.p")
ENCSR312AOC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR312AOC.p")
ENCSR165FLF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR165FLF.p")
ENCSR111SCI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR111SCI.p")
ENCSR662TPD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR662TPD.p")
ENCSR998YMP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR998YMP.p")
ENCSR776XMP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR776XMP.p")
ENCSR396XOO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR396XOO.p")
ENCSR711LRI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR711LRI.p")
ENCSR826EHD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR826EHD.p")
ENCSR108ZIJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR108ZIJ.p")
ENCSR097PDV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR097PDV.p")
ENCSR693KAV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR693KAV.p")
ENCSR899EKD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR899EKD.p")
ENCSR526LFP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR526LFP.p")
ENCSR627TLB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR627TLB.p")
ENCSR472OOH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR472OOH.p")
ENCSR049KDI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR049KDI.p")
ENCSR551OZP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR551OZP.p")
ENCSR358SBF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR358SBF.p")
ENCSR711CNO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR711CNO.p")
ENCSR291RIX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR291RIX.p")
ENCSR148VQR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR148VQR.p")
ENCSR196TWX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR196TWX.p")
ENCSR575TNG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR575TNG.p")
ENCSR950SIX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR950SIX.p")
ENCSR452BKB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR452BKB.p")
ENCSR035ZCY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR035ZCY.p")
ENCSR982GYD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR982GYD.p")
ENCSR732MDJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR732MDJ.p")
ENCSR422JWZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR422JWZ.p")
ENCSR763LIT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR763LIT.p")
ENCSR552XCF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR552XCF.p")
ENCSR093CLX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR093CLX.p")
ENCSR403MNL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR403MNL.p")
ENCSR165VMQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR165VMQ.p")
ENCSR263VNQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR263VNQ.p")
ENCSR402GQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR402GQZ.p")
ENCSR794RPU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR794RPU.p")
ENCSR454ZHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR454ZHR.p")
ENCSR216UTT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR216UTT.p")
ENCSR201SNQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR201SNQ.p")
ENCSR773CKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_hvg/ENCSR773CKJ.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_hvg_norm/Mus_musculus:adrenal_gland.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR402PHR,ENCSR157YJO,ENCSR021YYL,ENCSR988AMY,ENCSR908CQZ,ENCSR749GDE,ENCSR224OUG,ENCSR356VJZ,ENCSR312AOC,ENCSR165FLF,ENCSR111SCI,ENCSR662TPD,ENCSR998YMP,ENCSR776XMP,ENCSR396XOO,ENCSR711LRI,ENCSR826EHD,ENCSR108ZIJ,ENCSR097PDV,ENCSR693KAV,ENCSR899EKD,ENCSR526LFP,ENCSR627TLB,ENCSR472OOH,ENCSR049KDI,ENCSR551OZP,ENCSR358SBF,ENCSR711CNO,ENCSR291RIX,ENCSR148VQR,ENCSR196TWX,ENCSR575TNG,ENCSR950SIX,ENCSR452BKB,ENCSR035ZCY,ENCSR982GYD,ENCSR732MDJ,ENCSR422JWZ,ENCSR763LIT,ENCSR552XCF,ENCSR093CLX,ENCSR403MNL,ENCSR165VMQ,ENCSR263VNQ,ENCSR402GQZ,ENCSR794RPU,ENCSR454ZHR,ENCSR216UTT,ENCSR201SNQ,ENCSR773CKJ], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

