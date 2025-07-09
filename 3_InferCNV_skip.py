#%% Comprehensive InferCNV Analysis Pipeline for Epithelial Cells
import os
import gc
import scanpy as sc
import pandas as pd
import numpy as np
import infercnvpy as cnv
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

from infercnv_functions import *

print("🧬 Starting Comprehensive InferCNV Analysis Pipeline")
print("="*70)

# Set paths
workDir = "/mnt/public/lyx/CRC_dMMR"
figurePath = os.path.join(workDir, "Figures")
save_path = os.path.join(workDir, "Data", "processed", "scRNA", "integration_results")
infercnv_path = os.path.join(workDir, "Data", "processed", "scRNA", "inferCNV_results")

# Create inferCNV directory
if not os.path.exists(infercnv_path):
    os.makedirs(infercnv_path)

# Create subdirectories for different analyses
raw_cnv_path = os.path.join(infercnv_path, "raw_epithelial")
harmony_cnv_path = os.path.join(infercnv_path, "harmony_epithelial")
comparison_path = os.path.join(infercnv_path, "comparison")

for path in [raw_cnv_path, harmony_cnv_path, comparison_path]:
    if not os.path.exists(path):
        os.makedirs(path)

#%% Step 1: Subset Epithelial Cells from Integration Data
print("\n📋 Step 1: Subsetting Epithelial Cells from Integration Data")
# print("="*60)
print("="*80)

# Load integrated data if not already loaded
if 'final_integrated' not in locals():
    final_integrated = sc.read_h5ad(f"{save_path}/major_anno_all.h5ad")

# Check major cell types
print("Major cell type distribution:")
print(final_integrated.obs['Major_type'].value_counts())

# Subset epithelial cells
epithelial_mask = final_integrated.obs['Major_type'] == 'Epithelial'
# epithelial_cells = final_integrated[epithelial_mask].copy()
epithelial_cells = sc.AnnData(
        X=final_integrated[epithelial_mask].raw.X,
        obs=final_integrated[epithelial_mask].obs,
        var=final_integrated.raw.var
    )

print(f"\n✅ Extracted epithelial cells:")
print(f"   Total epithelial cells: {epithelial_cells.n_obs:,}")
print(f"   Features: {epithelial_cells.n_vars:,}")

# Check tissue type distribution in epithelial cells
print(f"\nTissue type distribution in epithelial cells:")
if 'Tissue_Type' in epithelial_cells.obs.columns:
    print(epithelial_cells.obs['Tissue_Type'].value_counts())

# Check treatment stage distribution
print(f"\nTreatment stage distribution in epithelial cells:")
if 'Treatment_Stage' in epithelial_cells.obs.columns:
    print(epithelial_cells.obs['Treatment_Stage'].value_counts())

# Check batch distribution
print(f"\nBatch distribution in epithelial cells:")
print(epithelial_cells.obs['batch'].value_counts())

#%% Step 2: Prepare Data Copies for Analysis
print("\n📋 Step 2: Preparing Data Copies for Analysis")
print("="*60)

# Copy 1: Raw epithelial data (no additional batch correction)
epithelial_raw = epithelial_cells.copy()
print("✅ Created raw epithelial copy")

# Copy 2: Harmony-corrected epithelial data
epithelial_harmony = epithelial_cells.copy()
print("✅ Created harmony epithelial copy")

# Prepare annotations for both datasets
epithelial_raw = prepare_cnv_annotations(epithelial_raw)
epithelial_harmony = prepare_cnv_annotations(epithelial_harmony)

#%% Step 3: Quality Control for InferCNV
print("\n📋 Step 3: Quality Control for InferCNV")
print("="*60)

# Apply QC to both datasets
epithelial_raw = infercnv_qc(epithelial_raw, max_mt_pct=75)
epithelial_harmony = infercnv_qc(epithelial_harmony, max_mt_pct=75)

#%% Step 4: InferCNV on Raw Data (No Additional Batch Correction)
print("\n🧬 Step 4: InferCNV Analysis on Raw Epithelial Data")
print("="*60)

# Run InferCNV on raw data
epithelial_raw_cnv = run_infercnv_analysis(
    adata = epithelial_raw, 
    output_path = raw_cnv_path, 
    analysis_name = "raw_epithelial",
    gtf_path="./GENCODE/HS.gencode.v42.annotation.gtf.gz"  # Specify your GTF path here if available
)

#%% Step 5: Additional Harmony Batch Correction for Epithelial Cells
print("\n🎯 Step 5: Additional Harmony Correction for Epithelial Cells")
print("="*60)

# Additional harmony correction specifically for epithelial cells
# (Note: Your data already has harmony correction, but this adds epithelial-specific correction)

# First, recalculate PCA for epithelial cells only
sc.pp.highly_variable_genes(epithelial_harmony, n_top_genes=2000)
sc.pp.scale(epithelial_harmony)
sc.tl.pca(epithelial_harmony, n_comps=50)

# Apply harmony integration specifically to epithelial cells
try:
    sc.external.pp.harmony_integrate(
        epithelial_harmony, 
        key='batch',
        basis='X_pca',
        adjusted_basis='X_pca_harmony_epi'
    )
    print("✅ Additional harmony correction completed for epithelial cells")
except Exception as e:
    print(f"⚠️  Harmony correction failed, using existing harmony: {str(e)}")
    # Use existing harmony coordinates
    epithelial_harmony.obsm['X_pca_harmony_epi'] = epithelial_harmony.obsm['X_pca_harmony']

#%% Step 6: InferCNV on Harmony-Corrected Data
print("\n🧬 Step 6: InferCNV Analysis on Harmony-Corrected Epithelial Data")
print("="*60)

# Run InferCNV on harmony-corrected data
epithelial_harmony_cnv = run_infercnv_analysis(
    epithelial_harmony, 
    harmony_cnv_path, 
    "harmony_epithelial",
    gtf_path="./GENCODE/HS.gencode.v42.annotation.gtf.gz"  # Specify your GTF path here if available
)

#%% Step 7: Save Results
print("\n💾 Step 7: Saving Results")
print("="*60)

# Save raw CNV results
epithelial_raw_cnv.write(os.path.join(raw_cnv_path, "epithelial_raw_cnv.h5ad"))
epithelial_raw_cnv.obs.to_csv(os.path.join(raw_cnv_path, "epithelial_raw_cnv_metadata.csv"))

# Save harmony CNV results  
epithelial_harmony_cnv.write(os.path.join(harmony_cnv_path, "epithelial_harmony_cnv.h5ad"))
epithelial_harmony_cnv.obs.to_csv(os.path.join(harmony_cnv_path, "epithelial_harmony_cnv_metadata.csv"))

print("✅ Results saved successfully")

#%% Step 8: Compare Results
print("\n🔍 Step 8: Comparing Raw vs Harmony InferCNV Results")
print("="*60)


# Perform comparison
comparison_results = compare_cnv_results(epithelial_raw_cnv, epithelial_harmony_cnv, comparison_path)

#%% Step 9: Final Summary and Recommendations
print("\n📊 Step 9: Final Summary and Recommendations")
print("="*60)

print(f"✅ InferCNV Analysis Completed Successfully!")
print(f"\n📁 Results saved to:")
print(f"   Raw analysis: {raw_cnv_path}")
print(f"   Harmony analysis: {harmony_cnv_path}")
print(f"   Comparison: {comparison_path}")

print(f"\n📋 Summary Statistics:")
print(f"   Raw epithelial cells analyzed: {epithelial_raw_cnv.n_obs:,}")
print(f"   Harmony epithelial cells analyzed: {epithelial_harmony_cnv.n_obs:,}")

if 'cnv_score' in epithelial_raw_cnv.obs.columns:
    print(f"   Raw CNV score range: {epithelial_raw_cnv.obs['cnv_score'].min():.3f} to {epithelial_raw_cnv.obs['cnv_score'].max():.3f}")

if 'cnv_score' in epithelial_harmony_cnv.obs.columns:
    print(f"   Harmony CNV score range: {epithelial_harmony_cnv.obs['cnv_score'].min():.3f} to {epithelial_harmony_cnv.obs['cnv_score'].max():.3f}")

print(f"\n🎯 Next Steps for Phase II:")
print(f"   1. Examine CNV scores to identify malignant cells")
print(f"   2. Choose raw vs harmony based on comparison results")
print(f"   3. Filter malignant epithelial cells for subtyping")
print(f"   4. Proceed with NMF clustering for Subtype X discovery")

print(f"\n🎉 InferCNV pipeline completed successfully!")