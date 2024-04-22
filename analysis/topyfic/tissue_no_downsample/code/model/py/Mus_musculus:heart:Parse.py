import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Mus_musculus:heart:Parse')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Mus_musculus:heart:Parse')

ENCSR729BRE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR729BRE.p")
ENCSR266RTS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR266RTS.p")
ENCSR442NQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR442NQZ.p")
ENCSR440PMT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR440PMT.p")
ENCSR938YKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR938YKJ.p")
ENCSR041KFB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR041KFB.p")
ENCSR906KCJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR906KCJ.p")
ENCSR764CYB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR764CYB.p")
ENCSR594LDC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR594LDC.p")
ENCSR860LQG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR860LQG.p")
ENCSR303EIR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR303EIR.p")
ENCSR617CNE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR617CNE.p")
ENCSR305SPL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR305SPL.p")
ENCSR147TEP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR147TEP.p")
ENCSR122DEZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR122DEZ.p")
ENCSR247PQY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR247PQY.p")
ENCSR868XSG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR868XSG.p")
ENCSR975YAQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR975YAQ.p")
ENCSR125HRM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR125HRM.p")
ENCSR277UKO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR277UKO.p")
ENCSR237OCP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR237OCP.p")
ENCSR688APT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR688APT.p")
ENCSR559AAG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR559AAG.p")
ENCSR337DKH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR337DKH.p")
ENCSR809WJL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR809WJL.p")
ENCSR573UKV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR573UKV.p")
ENCSR743OKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR743OKJ.p")
ENCSR455YDL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR455YDL.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_no_downsample/Mus_musculus:heart:Parse.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR729BRE,ENCSR266RTS,ENCSR442NQZ,ENCSR440PMT,ENCSR938YKJ,ENCSR041KFB,ENCSR906KCJ,ENCSR764CYB,ENCSR594LDC,ENCSR860LQG,ENCSR303EIR,ENCSR617CNE,ENCSR305SPL,ENCSR147TEP,ENCSR122DEZ,ENCSR247PQY,ENCSR868XSG,ENCSR975YAQ,ENCSR125HRM,ENCSR277UKO,ENCSR237OCP,ENCSR688APT,ENCSR559AAG,ENCSR337DKH,ENCSR809WJL,ENCSR573UKV,ENCSR743OKJ,ENCSR455YDL], data=data)

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


