import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random
import os

os.mkdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/res/Mus_musculus:gastrocnemius')
os.chdir('/home/whou10/scratch4/whou10/encode4/topyfic/tissue/res/Mus_musculus:gastrocnemius')

ENCSR478PUL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR478PUL.p")
ENCSR021MHR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR021MHR.p")
ENCSR438FRT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR438FRT.p")
ENCSR835ZXV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR835ZXV.p")
ENCSR188GDL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR188GDL.p")
ENCSR163GJZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR163GJZ.p")
ENCSR125XUV=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR125XUV.p")
ENCSR771MOE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR771MOE.p")
ENCSR627PLH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR627PLH.p")
ENCSR930FUM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR930FUM.p")
ENCSR554FGW=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR554FGW.p")
ENCSR667UMR=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR667UMR.p")
ENCSR355CTS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR355CTS.p")
ENCSR581NLU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR581NLU.p")
ENCSR627XNK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR627XNK.p")
ENCSR823CLM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR823CLM.p")
ENCSR488XVC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR488XVC.p")
ENCSR459ESG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR459ESG.p")
ENCSR576UHY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR576UHY.p")
ENCSR047KXM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR047KXM.p")
ENCSR927BSJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR927BSJ.p")
ENCSR019BZL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR019BZL.p")
ENCSR696NAL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR696NAL.p")
ENCSR328JOY=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR328JOY.p")
ENCSR253YLA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR253YLA.p")
ENCSR653HPZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR653HPZ.p")
ENCSR462NDO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR462NDO.p")
ENCSR498SXJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR498SXJ.p")
ENCSR399TQO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR399TQO.p")
ENCSR577OFZ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR577OFZ.p")
ENCSR966QJD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR966QJD.p")
ENCSR580UNE=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR580UNE.p")
ENCSR137UJM=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR137UJM.p")
ENCSR150KKS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR150KKS.p")
ENCSR279GBC=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR279GBC.p")
ENCSR480RWU=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR480RWU.p")
ENCSR928UPG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR928UPG.p")
ENCSR419OCT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR419OCT.p")
ENCSR921SBQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR921SBQ.p")
ENCSR650QSN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR650QSN.p")
ENCSR650VIQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR650VIQ.p")
ENCSR728NFJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR728NFJ.p")
ENCSR161BNJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR161BNJ.p")
ENCSR514XCG=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR514XCG.p")
ENCSR961HBK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR961HBK.p")
ENCSR305VTL=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR305VTL.p")
ENCSR908PNN=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR908PNN.p")
ENCSR496OYJ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR496OYJ.p")
ENCSR974ODT=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR974ODT.p")
ENCSR187HKH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR187HKH.p")
ENCSR644SME=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR644SME.p")
ENCSR790DPF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR790DPF.p")
ENCSR138QIF=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR138QIF.p")
ENCSR182SMX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR182SMX.p")
ENCSR094YHS=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR094YHS.p")
ENCSR160AKB=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR160AKB.p")
ENCSR800QGK=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR800QGK.p")
ENCSR358NKP=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR358NKP.p")
ENCSR548JRH=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR548JRH.p")
ENCSR231VJO=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR231VJO.p")
ENCSR816PAX=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR816PAX.p")
ENCSR026VDA=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR026VDA.p")
ENCSR083EPD=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR083EPD.p")
ENCSR425UJQ=Topyfic.read_train("/home/whou10/scratch4/whou10/encode4/topyfic/train/res/ENCSR425UJQ.p")

def normalization(mtx):
    pf = mtx.sum(axis=1).ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    pf = log1p_pf.sum(axis=1).ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    return pf_log1p_pf

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_tissue/Mus_musculus:gastrocnemius.h5ad')
data.X = normalization(data.X)
data.X = np.round_(data.X)
    
top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[ENCSR478PUL,ENCSR021MHR,ENCSR438FRT,ENCSR835ZXV,ENCSR188GDL,ENCSR163GJZ,ENCSR125XUV,ENCSR771MOE,ENCSR627PLH,ENCSR930FUM,ENCSR554FGW,ENCSR667UMR,ENCSR355CTS,ENCSR581NLU,ENCSR627XNK,ENCSR823CLM,ENCSR488XVC,ENCSR459ESG,ENCSR576UHY,ENCSR047KXM,ENCSR927BSJ,ENCSR019BZL,ENCSR696NAL,ENCSR328JOY,ENCSR253YLA,ENCSR653HPZ,ENCSR462NDO,ENCSR498SXJ,ENCSR399TQO,ENCSR577OFZ,ENCSR966QJD,ENCSR580UNE,ENCSR137UJM,ENCSR150KKS,ENCSR279GBC,ENCSR480RWU,ENCSR928UPG,ENCSR419OCT,ENCSR921SBQ,ENCSR650QSN,ENCSR650VIQ,ENCSR728NFJ,ENCSR161BNJ,ENCSR514XCG,ENCSR961HBK,ENCSR305VTL,ENCSR908PNN,ENCSR496OYJ,ENCSR974ODT,ENCSR187HKH,ENCSR644SME,ENCSR790DPF,ENCSR138QIF,ENCSR182SMX,ENCSR094YHS,ENCSR160AKB,ENCSR800QGK,ENCSR358NKP,ENCSR548JRH,ENCSR231VJO,ENCSR816PAX,ENCSR026VDA,ENCSR083EPD,ENCSR425UJQ], data=data)

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

