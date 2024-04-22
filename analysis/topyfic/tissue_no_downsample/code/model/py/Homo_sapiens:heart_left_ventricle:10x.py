import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Homo_sapiens:heart_left_ventricle:10x')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Homo_sapiens:heart_left_ventricle:10x')

ENCSR994VEG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR994VEG.p")
ENCSR067BOK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR067BOK.p")
ENCSR176WWW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR176WWW.p")
ENCSR398YBK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR398YBK.p")
ENCSR980OCK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR980OCK.p")
ENCSR405YKM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR405YKM.p")
ENCSR481QQR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR481QQR.p")
ENCSR755CIF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR755CIF.p")
ENCSR727OYO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR727OYO.p")
ENCSR906MRL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR906MRL.p")
ENCSR084XKX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR084XKX.p")
ENCSR345CVL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR345CVL.p")
ENCSR349AHE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR349AHE.p")
ENCSR962JKS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR962JKS.p")
ENCSR681GPY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR681GPY.p")
ENCSR273JWD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR273JWD.p")
ENCSR630LZS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR630LZS.p")
ENCSR654MFX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR654MFX.p")
ENCSR157FDD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR157FDD.p")
ENCSR814LMX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR814LMX.p")
ENCSR231FNL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR231FNL.p")
ENCSR540DHJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR540DHJ.p")
ENCSR488UUT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR488UUT.p")
ENCSR008CVR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR008CVR.p")
ENCSR056QLB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR056QLB.p")
ENCSR489URW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR489URW.p")
ENCSR204RHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR204RHR.p")
ENCSR328GTN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR328GTN.p")
ENCSR788SNY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR788SNY.p")
ENCSR991LHO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR991LHO.p")
ENCSR899GYX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR899GYX.p")
ENCSR203YOV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR203YOV.p")
ENCSR012APQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR012APQ.p")
ENCSR259VOY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR259VOY.p")
ENCSR455MGH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR455MGH.p")
ENCSR076ZLE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR076ZLE.p")
ENCSR763BII=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR763BII.p")
ENCSR190TRK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR190TRK.p")
ENCSR002SMQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR002SMQ.p")
ENCSR485GOL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR485GOL.p")
ENCSR762LML=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR762LML.p")
ENCSR801DHT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR801DHT.p")
ENCSR439ZVQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR439ZVQ.p")
ENCSR352DXB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR352DXB.p")
ENCSR751BHQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR751BHQ.p")
ENCSR919ENI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR919ENI.p")
ENCSR138JCM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR138JCM.p")
ENCSR175TRJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR175TRJ.p")
ENCSR753YOZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR753YOZ.p")
ENCSR237HWJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR237HWJ.p")
ENCSR093GXF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR093GXF.p")
ENCSR777RUZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR777RUZ.p")
ENCSR851JBE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR851JBE.p")
ENCSR085XEW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR085XEW.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_no_downsample/Homo_sapiens:heart_left_ventricle:10x.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR994VEG,ENCSR067BOK,ENCSR176WWW,ENCSR398YBK,ENCSR980OCK,ENCSR405YKM,ENCSR481QQR,ENCSR755CIF,ENCSR727OYO,ENCSR906MRL,ENCSR084XKX,ENCSR345CVL,ENCSR349AHE,ENCSR962JKS,ENCSR681GPY,ENCSR273JWD,ENCSR630LZS,ENCSR654MFX,ENCSR157FDD,ENCSR814LMX,ENCSR231FNL,ENCSR540DHJ,ENCSR488UUT,ENCSR008CVR,ENCSR056QLB,ENCSR489URW,ENCSR204RHR,ENCSR328GTN,ENCSR788SNY,ENCSR991LHO,ENCSR899GYX,ENCSR203YOV,ENCSR012APQ,ENCSR259VOY,ENCSR455MGH,ENCSR076ZLE,ENCSR763BII,ENCSR190TRK,ENCSR002SMQ,ENCSR485GOL,ENCSR762LML,ENCSR801DHT,ENCSR439ZVQ,ENCSR352DXB,ENCSR751BHQ,ENCSR919ENI,ENCSR138JCM,ENCSR175TRJ,ENCSR753YOZ,ENCSR237HWJ,ENCSR093GXF,ENCSR777RUZ,ENCSR851JBE,ENCSR085XEW], data=data)

top_model.save_topModel()
analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)
analysis_top_model.save_analysis()

metadata = ['celltype']
analysis_top_model.TopicTraitRelationshipHeatmap(metaData=metadata,
                                                 save=True,
                                                 show=True,
                                                 file_name='figures/topic-traitRelationships')

label = dict()
for i in range(analysis_top_model.top_model.N):
    key = f"Topic_{i+1}"
    value = f"sc _{i+1}"
    label[key] = value

analysis_top_model.average_cell_participation(label=label,
                                              figsize=(15,5),
                                              file_name="figures/average_cell_participation")

file_name = f"figures/pieChart_dist_topics_subtypes_RNA"
analysis_top_model.pie_structure_Chart(level='celltype',
                                       save=True,
                                       show=True,
                                       file_name=file_name)
counts=data.obs.value_counts()
counts=counts[counts > 50]
ct=counts.index.tolist()
ct = [item[0] for item in ct]

file_name = f"figures/structurePlot_dist_topics_cellType_RNA_cluster_hierarchical_State"
analysis_top_model.structure_plot(level='celltype',
                                  category=ct, 
                                  save=True,
                                  show=True,
                                  figsize=(100,20),
                                  file_name=file_name)

# gene_weights = top_model.get_gene_weights()
# gene_weights
# 
# GSEA_topic1 = top_model.topics['Topic_1'].GSEA(file_name=f"figures/GSEA_topic_1",
#                                                  verbose=True)


