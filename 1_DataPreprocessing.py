import os
import gc

import scanpy as sc

import pandas as pd
import numpy as np

## Function to load and process scRNA-seq datasets
from preprocessing_functions import *

#%% Main here
## Set workDir
workDir = "/mnt/public/lyx/CRC_dMMR"
figurePath = os.path.join(workDir, "Figures")

Raw_scrna_data_Path = os.path.join(workDir, "Data","raw", "scRNA")
Processed_scrna_data_Path = os.path.join(workDir, "Data", "processed", "scRNA")

if not os.path.exists(figurePath):
    os.makedirs(figurePath)

if not os.path.exists(Processed_scrna_data_Path):
    os.makedirs(Processed_scrna_data_Path)

## View all data
dataFiles = os.listdir(Raw_scrna_data_Path)

preprocessed_datasets = {}

## Preprocess each dataset
for study_name in dataFiles:
    # print(f"Study: {study_name} with {os.listdir(os.path.join(Raw_scrna_data_Path, study_name))} files")

    ## Load and process each dataset
    data_path = os.path.join(Raw_scrna_data_Path, study_name)
    clinical_data_path = os.path.join(data_path,"clinical_data.csv")

    ## Check the Processed data existing
    if os.path.exists(os.path.join(Processed_scrna_data_Path,f"preprocessed_{study_name}.h5ad")):
        print(f"Preprocessed data for {study_name} already exists. Skipping preprocessing.")
        
        ## Read and combine preprocessed data
        preprocessed_adata = sc.read(os.path.join(Processed_scrna_data_Path,f"preprocessed_{study_name}.h5ad"))
        preprocessed_datasets[study_name] = preprocessed_adata

        del preprocessed_adata
        gc.collect()

    else:
        ## Check the raw file existing
        if os.path.exists(os.path.join(Processed_scrna_data_Path,f"raw_{study_name}.h5ad")):
            print(f"Loading dataset: {study_name} from {data_path}")
            adata = sc.read(os.path.join(Processed_scrna_data_Path,f"raw_{study_name}.h5ad"))
            print(f"✅ Successfully loaded {study_name} from processed data")

        else:
            print(f"Raw data for {study_name} not exist. Start loading step.")
            adata = load_scrna_dataset(study_name, data_path)

            if adata is not None:
                print(f"✅ Successfully loaded {study_name}")
            else:
                print(f"❌ Failed to load {study_name}")

            ## Save raw data
            adata.write(os.path.join(Processed_scrna_data_Path,f"raw_{study_name}.h5ad"))

        ## Start Preprocessing
        try:
            preprocessed_adata = preprocess_single_dataset(
                adata.copy(), 
                study_name, 
                clinical_data_path
            )

            ## Combine preprocessed data
            preprocessed_datasets[study_name] = preprocessed_adata
            
            # Print summary
            print(f"\n📊 {study_name} Summary:")
            print(f"  Final shape: {preprocessed_adata.shape}")
            print(f"  Doublets flagged: {preprocessed_adata.obs['predicted_doublet'].sum()}")
            print(f"  Mean genes/cell: {preprocessed_adata.obs['n_genes_by_counts'].mean():.0f}")
            print(f"  Mean MT%: {preprocessed_adata.obs['pct_counts_mt'].mean():.2f}")
            
        except Exception as e:
            print(f"❌ Failed to preprocess {study_name}: {str(e)}")
            continue
    
        ## Save data
        preprocessed_adata.write(os.path.join(Processed_scrna_data_Path,f"preprocessed_{study_name}.h5ad"))
        del adata, preprocessed_adata
        gc.collect()

## Perform integration
#%%
merged_adata = merge_datasets_for_integration(preprocessed_datasets)
n_top_genes=2000

# Phase 1: Quick method comparison (100K cells)
quick_results, quick_metrics, quick_fig = quick_method_comparison(
    merged_adata, 
    batch_key='batch',
    n_sample=80000,  # Fast comparison
    n_top_genes=2000
)

print("🔍 Quick comparison metrics:")
print(quick_metrics)

print(quick_fig)
quick_fig.savefig(
    os.path.join(figurePath, "quick_method_comparison.pdf"), 
    dpi=300, bbox_inches='tight')

# Phase 2: Integration without downsampling
print(f"\n🎯 PHASE 2: Detailed Integration")
print("="*60)

# Prepare data with larger sample
final_integrated = prepare_for_integration(
    merged_adata, 
    n_top_genes=n_top_genes,
    downsample_n=None,
    batch_key='batch'
)

sc.external.pp.harmony_integrate(
    final_integrated, 
    key='batch',
    basis='X_pca',
    adjusted_basis='X_pca_harmony'
)

# Save final result
save_path = os.path.join(workDir, "Data", "processed", "scRNA", "integration_results")
if not os.path.exists(save_path):
    os.makedirs(save_path)

final_integrated.write(f"{save_path}/final_integrated.h5ad")