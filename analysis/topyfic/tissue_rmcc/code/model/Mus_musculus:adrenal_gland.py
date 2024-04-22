import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Mus_musculus:adrenal_gland')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_rmcc/res/Mus_musculus:adrenal_gland')

ENCSR402PHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR402PHR.p")
ENCSR157YJO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR157YJO.p")
ENCSR021YYL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR021YYL.p")
ENCSR988AMY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR988AMY.p")
ENCSR908CQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR908CQZ.p")
ENCSR749GDE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR749GDE.p")
ENCSR224OUG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR224OUG.p")
ENCSR356VJZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR356VJZ.p")
ENCSR312AOC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR312AOC.p")
ENCSR165FLF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR165FLF.p")
ENCSR111SCI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR111SCI.p")
ENCSR662TPD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR662TPD.p")
ENCSR998YMP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR998YMP.p")
ENCSR776XMP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR776XMP.p")
ENCSR396XOO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR396XOO.p")
ENCSR711LRI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR711LRI.p")
ENCSR826EHD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR826EHD.p")
ENCSR108ZIJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR108ZIJ.p")
ENCSR097PDV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR097PDV.p")
ENCSR693KAV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR693KAV.p")
ENCSR899EKD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR899EKD.p")
ENCSR526LFP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR526LFP.p")
ENCSR627TLB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR627TLB.p")
ENCSR472OOH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR472OOH.p")
ENCSR049KDI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR049KDI.p")
ENCSR551OZP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR551OZP.p")
ENCSR358SBF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR358SBF.p")
ENCSR711CNO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR711CNO.p")
ENCSR291RIX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR291RIX.p")
ENCSR148VQR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR148VQR.p")
ENCSR196TWX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR196TWX.p")
ENCSR575TNG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR575TNG.p")
ENCSR950SIX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR950SIX.p")
ENCSR452BKB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR452BKB.p")
ENCSR035ZCY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR035ZCY.p")
ENCSR982GYD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR982GYD.p")
ENCSR732MDJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR732MDJ.p")
ENCSR422JWZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR422JWZ.p")
ENCSR763LIT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR763LIT.p")
ENCSR552XCF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR552XCF.p")
ENCSR093CLX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR093CLX.p")
ENCSR403MNL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR403MNL.p")
ENCSR165VMQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR165VMQ.p")
ENCSR263VNQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR263VNQ.p")
ENCSR402GQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR402GQZ.p")
ENCSR794RPU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR794RPU.p")
ENCSR454ZHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR454ZHR.p")
ENCSR216UTT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR216UTT.p")
ENCSR201SNQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR201SNQ.p")
ENCSR773CKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/ENCSR773CKJ.p")

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_rmcc_norm/Mus_musculus:adrenal_gland.h5ad')
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR402PHR,ENCSR157YJO,ENCSR021YYL,ENCSR988AMY,ENCSR908CQZ,ENCSR749GDE,ENCSR224OUG,ENCSR356VJZ,ENCSR312AOC,ENCSR165FLF,ENCSR111SCI,ENCSR662TPD,ENCSR998YMP,ENCSR776XMP,ENCSR396XOO,ENCSR711LRI,ENCSR826EHD,ENCSR108ZIJ,ENCSR097PDV,ENCSR693KAV,ENCSR899EKD,ENCSR526LFP,ENCSR627TLB,ENCSR472OOH,ENCSR049KDI,ENCSR551OZP,ENCSR358SBF,ENCSR711CNO,ENCSR291RIX,ENCSR148VQR,ENCSR196TWX,ENCSR575TNG,ENCSR950SIX,ENCSR452BKB,ENCSR035ZCY,ENCSR982GYD,ENCSR732MDJ,ENCSR422JWZ,ENCSR763LIT,ENCSR552XCF,ENCSR093CLX,ENCSR403MNL,ENCSR165VMQ,ENCSR263VNQ,ENCSR402GQZ,ENCSR794RPU,ENCSR454ZHR,ENCSR216UTT,ENCSR201SNQ,ENCSR773CKJ], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

