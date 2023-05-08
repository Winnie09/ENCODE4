import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random


def normalization(mtx):
    pf = mtx.sum(axis=1).A.ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    
    pf = log1p_pf.sum(axis=1).A.ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    
    return pf_log1p_pf


data = sc.read('/home/zji4/data-hji7/zji/encode/data/scrna/anndata/h5ad/ENCSR966DDY.rds.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)

train = Topyfic.Train(name="model",k=15,n_runs=100)
train.run_LDA_models(data, n_jobs=5, n_thread=5)

train.save_train()
train = Topyfic.read_train("train_model.p")


top_model, clustering = Topyfic.calculate_leiden_clustering(trains=[train], data=data)

analysis_top_model = Topyfic.Analysis(Top_model=top_model)

analysis_top_model.calculate_cell_participation(data=data)

metadata = ['ct']
analysis_top_model.TopicTraitRelationshipHeatmap(metaData=metadata,
                                                 save=True,
                                                 show=True,
                                                 file_name='figures/topic-traitRelationships')
                                                 
                                                 
cellType = ['alveolar type 2 cell','alveolar fibroblast','alveolar type 1 cell','macrophage','T cell']
file_name = "figures/pieChart_dist_topics_subtypes_RNA"
analysis_top_model.pie_structure_Chart(level='ct',
                                       category=cellType,
                                       save=True,
                                       show=True,
                                       file_name=file_name)


label = dict()
for i in range(analysis_top_model.top_model.N):
    key = f"Topic_{i+1}"
    value = f"sc _{i+1}"
    label[key] = value

analysis_top_model.average_cell_participation(label=label,
                                              figsize=(15,5),
                                              file_name="figures/average_cell_participation")
                                              
                                              
                                              
                                              

file_name = "figures/structurePlot_dist_topics_cellType_RNA_cluster_hierarchical_seurat_clusters"
analysis_top_model.structure_plot(level='ct',
                                  category=cellType, 
                                  order_cells=['ct'],
                                  save=True,
                                  show=True,
                                  figsize=(60,20),
                                  file_name=file_name)

gene_rank_weights = top_model.get_ranked_gene_weight()
gene_rank_weights.to_csv('generank.csv')


GSEA_topic1 = top_model.topics['Topic_1'].GSEA(file_name=f"figures/GSEA_topic_1",verbose=True)
GSEA_topic1.to_csv('GSEA_topic1.csv')
