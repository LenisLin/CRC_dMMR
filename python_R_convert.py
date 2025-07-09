import os
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse

## Set Path 
resultPath = "/home/lenislin/Experiment/projects/HorizontalProject/zsz_CRC/Results/anno"
dataPath = "/home/lenislin/Experiment/projects/HorizontalProject/zsz_CRC/Data"

## Load data
adata = sc.read_h5ad(os.path.join(dataPath,"majorAnno_adata.h5ad"))
adata_tumor = adata[adata.obs['major_celltype'] == 'Epithelial']

# Save the expression profile (X) as a sparse matrix in .mtx format
matrix = scipy.sparse.csr_matrix(np.array(adata_tumor.X))

dataPath = "/home/lenislin/Experiment/projects/HorizontalProject/zsz_CRC/Data/Step4"
scipy.io.mmwrite(os.path.join(dataPath,"expression_profile.mtx"), matrix)

# Save the metadata (obs) as a CSV file
adata_tumor.obs.to_csv(os.path.join(dataPath,'metadata.csv'))

# Save row names (var) and column names (obs) if needed
adata_tumor.var_names.to_series().to_csv(os.path.join(dataPath,'row_names.csv'), header=False)
adata_tumor.obs_names.to_series().to_csv(os.path.join(dataPath,'column_names.csv'), header=False)

