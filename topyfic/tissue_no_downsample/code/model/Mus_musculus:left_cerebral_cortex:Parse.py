import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Mus_musculus:left_cerebral_cortex:Parse')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Mus_musculus:left_cerebral_cortex:Parse')

ENCSR490RCS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR490RCS.p")
ENCSR037LFG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR037LFG.p")
ENCSR827AOS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR827AOS.p")
ENCSR527LIA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR527LIA.p")
ENCSR084RQP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR084RQP.p")
ENCSR118FCJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR118FCJ.p")
ENCSR611QJH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR611QJH.p")
ENCSR248MFG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR248MFG.p")
ENCSR450IJU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR450IJU.p")
ENCSR254SLA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR254SLA.p")
ENCSR839BHY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR839BHY.p")
ENCSR793HMA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR793HMA.p")
ENCSR308DHB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR308DHB.p")
ENCSR883XZK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR883XZK.p")
ENCSR632SJM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR632SJM.p")
ENCSR220GBV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR220GBV.p")
ENCSR257BOU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR257BOU.p")
ENCSR386PPF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR386PPF.p")
ENCSR081SXM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR081SXM.p")
ENCSR659DFF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR659DFF.p")
ENCSR964OQW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR964OQW.p")
ENCSR544MGC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR544MGC.p")
ENCSR302VFA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR302VFA.p")
ENCSR132KMR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR132KMR.p")
ENCSR178PKG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR178PKG.p")
ENCSR250SCW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR250SCW.p")
ENCSR092RYV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR092RYV.p")
ENCSR134FSZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR134FSZ.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_no_downsample/Mus_musculus:left_cerebral_cortex:Parse.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR490RCS,ENCSR037LFG,ENCSR827AOS,ENCSR527LIA,ENCSR084RQP,ENCSR118FCJ,ENCSR611QJH,ENCSR248MFG,ENCSR450IJU,ENCSR254SLA,ENCSR839BHY,ENCSR793HMA,ENCSR308DHB,ENCSR883XZK,ENCSR632SJM,ENCSR220GBV,ENCSR257BOU,ENCSR386PPF,ENCSR081SXM,ENCSR659DFF,ENCSR964OQW,ENCSR544MGC,ENCSR302VFA,ENCSR132KMR,ENCSR178PKG,ENCSR250SCW,ENCSR092RYV,ENCSR134FSZ], data=data)

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
