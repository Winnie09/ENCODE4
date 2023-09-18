import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/res/Mus_musculus:layer_of_hippocampus')
os.chdir('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/res/Mus_musculus:layer_of_hippocampus')

ENCSR944XPC=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR944XPC.p")
ENCSR995ZKH=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR995ZKH.p")
ENCSR539XHZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR539XHZ.p")
ENCSR987SOK=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR987SOK.p")
ENCSR835JFN=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR835JFN.p")
ENCSR968RPR=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR968RPR.p")
ENCSR742NWO=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR742NWO.p")
ENCSR075HVN=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR075HVN.p")
ENCSR079GGS=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR079GGS.p")
ENCSR846QIZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR846QIZ.p")
ENCSR483SSM=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR483SSM.p")
ENCSR642QUS=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR642QUS.p")
ENCSR569LJC=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR569LJC.p")
ENCSR399EUZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR399EUZ.p")
ENCSR967YET=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR967YET.p")
ENCSR940VCU=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR940VCU.p")
ENCSR203XJO=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR203XJO.p")
ENCSR435CMT=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR435CMT.p")
ENCSR545URU=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR545URU.p")
ENCSR257VQO=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR257VQO.p")
ENCSR127LUZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR127LUZ.p")
ENCSR859XOI=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR859XOI.p")
ENCSR736JKX=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR736JKX.p")
ENCSR499ZYZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR499ZYZ.p")
ENCSR983TAK=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR983TAK.p")
ENCSR613VKV=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR613VKV.p")
ENCSR767FKE=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR767FKE.p")
ENCSR824CFI=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR824CFI.p")
ENCSR329JYG=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR329JYG.p")
ENCSR079CKL=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR079CKL.p")
ENCSR888TQC=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR888TQC.p")
ENCSR979MRY=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR979MRY.p")
ENCSR973WCI=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR973WCI.p")
ENCSR247QPJ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR247QPJ.p")
ENCSR730JFZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR730JFZ.p")
ENCSR167HUO=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR167HUO.p")
ENCSR869HKZ=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR869HKZ.p")
ENCSR926QCE=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR926QCE.p")
ENCSR692LZA=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR692LZA.p")
ENCSR320EFC=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR320EFC.p")
ENCSR232VMC=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR232VMC.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/zji4/data-hji7/zji/encode/data/scrna/anndata/h5ad_tissue/Mus_musculus:layer_of_hippocampus.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR944XPC,ENCSR995ZKH,ENCSR539XHZ,ENCSR987SOK,ENCSR835JFN,ENCSR968RPR,ENCSR742NWO,ENCSR075HVN,ENCSR079GGS,ENCSR846QIZ,ENCSR483SSM,ENCSR642QUS,ENCSR569LJC,ENCSR399EUZ,ENCSR967YET,ENCSR940VCU,ENCSR203XJO,ENCSR435CMT,ENCSR545URU,ENCSR257VQO,ENCSR127LUZ,ENCSR859XOI,ENCSR736JKX,ENCSR499ZYZ,ENCSR983TAK,ENCSR613VKV,ENCSR767FKE,ENCSR824CFI,ENCSR329JYG,ENCSR079CKL,ENCSR888TQC,ENCSR979MRY,ENCSR973WCI,ENCSR247QPJ,ENCSR730JFZ,ENCSR167HUO,ENCSR869HKZ,ENCSR926QCE,ENCSR692LZA,ENCSR320EFC,ENCSR232VMC], data=data)

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
