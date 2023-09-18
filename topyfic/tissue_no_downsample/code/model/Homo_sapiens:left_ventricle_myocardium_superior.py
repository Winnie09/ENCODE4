import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/res/Homo_sapiens:left_ventricle_myocardium_superior')
os.chdir('/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/tissue/res/Homo_sapiens:left_ventricle_myocardium_superior')

ENCSR079LVG=Topyfic.read_train("/home/zji4/data-hji7/zji/encode/data/scrna/topyfic/train/res/ENCSR079LVG.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/zji4/data-hji7/zji/encode/data/scrna/anndata/h5ad_tissue/Homo_sapiens:left_ventricle_myocardium_superior.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR079LVG], data=data)

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
