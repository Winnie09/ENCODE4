import Topyfic
import sys
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as sp
import random

data = sc.read('/home/whou10/scratch4/whou10/encode4/data/scrna/anndata/h5ad_sample_rmcc_norm/'+sys.argv[1]+'.h5ad')
data.X = np.round_(data.X)

train = Topyfic.Train(name="model",k=15,n_runs=100)
train.run_LDA_models(data, n_jobs=5, n_thread=5)

train.save_train('/home/whou10/scratch4/whou10/encode4/topyfic/train_rmcc/'+sys.argv[1])


