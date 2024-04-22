import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os





os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/res/Homo_sapiens:lung')

top_model = Topyfic.read_topModel("topModel_model.p")

analysis_top_model = Topyfic.read_analysis("analysis_model.p")


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
analysis_top_model.cell_participation.to_df().to_csv('figures/cell_participation.csv')
counts=analysis_top_model.cell_participation.obs.value_counts()
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

gene_weights = top_model.get_gene_weights()
gene_weights.to_csv('figures/gene_weights.csv')

gene_rank_weights = top_model.get_ranked_gene_weight()
gene_rank_weights.to_csv('figures/gene_rank_weights.csv')

for i in range(analysis_top_model.top_model.N):
    GSEA_topic = top_model.topics["Topic_" + str(i+1)].GSEA(file_name=f"figures/GSEA_topic_{i+1}")
    GSEA_topic.to_csv()


