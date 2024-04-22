import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Mus_musculus:adrenal_gland:Parse')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue_no_downsample/res/Mus_musculus:adrenal_gland:Parse')

ENCSR312AOC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR312AOC.p")
ENCSR165FLF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR165FLF.p")
ENCSR111SCI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR111SCI.p")
ENCSR662TPD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR662TPD.p")
ENCSR998YMP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR998YMP.p")
ENCSR776XMP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR776XMP.p")
ENCSR396XOO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR396XOO.p")
ENCSR711LRI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR711LRI.p")
ENCSR826EHD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR826EHD.p")
ENCSR108ZIJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR108ZIJ.p")
ENCSR097PDV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR097PDV.p")
ENCSR693KAV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR693KAV.p")
ENCSR899EKD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR899EKD.p")
ENCSR526LFP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR526LFP.p")
ENCSR627TLB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR627TLB.p")
ENCSR472OOH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR472OOH.p")
ENCSR049KDI=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR049KDI.p")
ENCSR551OZP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR551OZP.p")
ENCSR358SBF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR358SBF.p")
ENCSR711CNO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR711CNO.p")
ENCSR291RIX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR291RIX.p")
ENCSR148VQR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR148VQR.p")
ENCSR196TWX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR196TWX.p")
ENCSR575TNG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR575TNG.p")
ENCSR950SIX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR950SIX.p")
ENCSR452BKB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR452BKB.p")
ENCSR035ZCY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR035ZCY.p")
ENCSR982GYD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR982GYD.p")
ENCSR732MDJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR732MDJ.p")
ENCSR422JWZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR422JWZ.p")
ENCSR763LIT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR763LIT.p")
ENCSR552XCF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR552XCF.p")
ENCSR093CLX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR093CLX.p")
ENCSR403MNL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR403MNL.p")
ENCSR165VMQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR165VMQ.p")
ENCSR263VNQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR263VNQ.p")
ENCSR402GQZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR402GQZ.p")
ENCSR794RPU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR794RPU.p")
ENCSR454ZHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR454ZHR.p")
ENCSR216UTT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR216UTT.p")
ENCSR201SNQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR201SNQ.p")
ENCSR773CKJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train_no_downsample/res/ENCSR773CKJ.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue_norm/Mus_musculus:adrenal_gland:Parse.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR312AOC,ENCSR165FLF,ENCSR111SCI,ENCSR662TPD,ENCSR998YMP,ENCSR776XMP,ENCSR396XOO,ENCSR711LRI,ENCSR826EHD,ENCSR108ZIJ,ENCSR097PDV,ENCSR693KAV,ENCSR899EKD,ENCSR526LFP,ENCSR627TLB,ENCSR472OOH,ENCSR049KDI,ENCSR551OZP,ENCSR358SBF,ENCSR711CNO,ENCSR291RIX,ENCSR148VQR,ENCSR196TWX,ENCSR575TNG,ENCSR950SIX,ENCSR452BKB,ENCSR035ZCY,ENCSR982GYD,ENCSR732MDJ,ENCSR422JWZ,ENCSR763LIT,ENCSR552XCF,ENCSR093CLX,ENCSR403MNL,ENCSR165VMQ,ENCSR263VNQ,ENCSR402GQZ,ENCSR794RPU,ENCSR454ZHR,ENCSR216UTT,ENCSR201SNQ,ENCSR773CKJ], data=data)

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



