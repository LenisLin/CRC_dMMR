#%% Comprehensive Malignant Cells Analysis Pipeline
import os
import gc
import pickle
import scanpy as sc
import gseapy as gp
from gseapy import Msigdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

from sklearn.decomposition import NMF
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
import scipy.sparse

import warnings
warnings.filterwarnings('ignore')

from malignant_functions import *

print("🧬 Starting Comprehensive Malignant Cells Analysis Pipeline")
print("="*70)

# Set paths
workDir = "/mnt/public/lyx/CRC_dMMR"
save_path = os.path.join(workDir, "Results", "integration_results")
result_path = os.path.join(workDir, "Results","NMF_results")
figurePath = os.path.join(workDir, "Figures","NMF_results")

# create directories if they don't exist
os.makedirs(result_path, exist_ok=True)
os.makedirs(figurePath, exist_ok=True)

#%% Task 1: Extract Malignant-Only Cells with Original Expression
print("\n🎯 Task 1: Extracting Malignant-Only Cells")
print("="*60)

# Load integrated data if not already loaded
if 'final_integrated' not in locals():
    final_integrated = sc.read_h5ad(f"{save_path}/major_anno_all.h5ad")

# Check major cell types
print("Major cell type distribution:")
print(final_integrated.obs['Major_type'].value_counts())
print(final_integrated.obs['Sample_ID'].value_counts())

# Filter for tumor-origin epithelial cells only (malignant)
malignant_mask = (final_integrated.obs['Major_type'] == 'Epithelial') & \
                 (final_integrated.obs['Tissue_Type'] == 'tumor')

print(f"Total epithelial cells: {(final_integrated.obs['Major_type'] == 'Epithelial').sum()}")
print(f"Tumor epithelial cells: {malignant_mask.sum()}")

# Create malignant-only object with original expression matrix
malignant_cells = sc.AnnData(
    X=final_integrated[malignant_mask].raw.X,
    obs=final_integrated[malignant_mask].obs.copy(),
    var=final_integrated.raw.var.copy()
)
malignant_cells.obs['Major_type'] = 'Malignant'

print(f"\n✅ Extracted malignant cells:")
print(f"   Total malignant cells: {malignant_cells.n_obs:,}")
print(f"   Total features: {malignant_cells.n_vars:,}")
print(f"   Batch distribution:")
print(malignant_cells.obs['study'].value_counts())

# Save malignant-only object
malignant_cells.write(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))
print(f"✅ Saved malignant cells to: {result_path}/malignant_epithelial_cells.h5ad")

#%% Task 2: Downsample and Visualize Batch Effects
print("\n🎨 Task 2: Downsampling and Batch Visualization")
print("="*60)

if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

# Downsample for visualization (stratified by batch)
target_cells = min(15000, malignant_cells.n_obs)
if malignant_cells.n_obs > target_cells:
    # Stratified sampling by batch
    batch_counts = malignant_cells.obs['batch'].value_counts()
    batch_proportions = batch_counts / batch_counts.sum()
    target_per_batch = (batch_proportions * target_cells).round().astype(int)
    
    sampled_indices = []
    for batch in malignant_cells.obs['batch'].unique():
        batch_mask = malignant_cells.obs['batch'] == batch
        batch_indices = np.where(batch_mask)[0]
        
        n_target = min(target_per_batch[batch], len(batch_indices))
        sampled_batch_indices = np.random.choice(batch_indices, size=n_target, replace=False)
        sampled_indices.extend(sampled_batch_indices)
    
    malignant_subset = malignant_cells[sampled_indices].copy()
else:
    malignant_subset = malignant_cells.copy()

print(f"Downsampled to {malignant_subset.n_obs} cells for visualization")

# Basic preprocessing for visualization
sc.pp.normalize_total(malignant_subset, target_sum=1e4)
sc.pp.log1p(malignant_subset)
sc.pp.highly_variable_genes(malignant_subset, n_top_genes=2000)
sc.pp.scale(malignant_subset)
sc.tl.pca(malignant_subset)
sc.pp.neighbors(malignant_subset, n_neighbors=15)
sc.tl.umap(malignant_subset)

# Create batch visualization
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Batch effect
sc.pl.umap(malignant_subset, color='batch', ax=axes[0,0], show=False, 
          legend_loc='right margin', title='Batch Effect')

# Treatment Strategy
sc.pl.umap(malignant_subset, color='Treatment_Strategy', ax=axes[0,1], show=False,
          legend_loc='right margin', title='Treatment Strategy')

# Treatment Stage
if 'Treatment_Stage' in malignant_subset.obs.columns:
    sc.pl.umap(malignant_subset, color='Treatment_Stage', ax=axes[0,2], show=False,
              legend_loc='right margin', title='Treatment Stage')

# Patient ID (to check patient mixing)
sc.pl.umap(malignant_subset, color='Sample_ID', ax=axes[1,0], show=False,
          legend_loc=None, title='Patient Distribution')

# Quality metrics
sc.pl.umap(malignant_subset, color='n_genes_by_counts', ax=axes[1,1], show=False,
          title='Gene Counts', color_map='viridis')

sc.pl.umap(malignant_subset, color='pct_counts_mt', ax=axes[1,2], show=False,
          title='MT Percentage', color_map='Reds')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, "malignant_batch_effects.pdf"), dpi=300, bbox_inches='tight')
plt.show()

print("✅ Batch visualization saved")

#%% Task 3: Extract Metabolism-Associated Genes using MSigDB API
print("\n🧬 Task 3: Extracting Metabolism-Associated Genes from MSigDB")
print("="*60)

malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

# Get available gene sets from MSigDB
metabolism_gene_sets = get_metabolism_pathways()
stress_gene_sets = get_stress_resistance_pathways()
oncogenic_gene_sets = get_oncogenic_pathways()
metabolism_gene_sets.update(stress_gene_sets)
metabolism_gene_sets.update(oncogenic_gene_sets)

print(f"\nTotal gene sets collected: {len(metabolism_gene_sets)}")

# Create pathway coverage report
pathway_coverage_report = []
for pathway_name, genes in metabolism_gene_sets.items():
    available_genes = [g for g in genes if g in malignant_cells.var_names]
    coverage = len(available_genes) / len(genes) * 100 if genes else 0
    
    pathway_coverage_report.append({
        'pathway': pathway_name,
        'total_genes': len(genes),
        'available_genes': len(available_genes),
        'coverage_percent': coverage,
        'available_gene_list': available_genes
    })

# Convert to DataFrame for easy viewing
coverage_df = pd.DataFrame(pathway_coverage_report)
coverage_df = coverage_df.sort_values('coverage_percent', ascending=False)

print(f"\nTop 10 best covered pathways:")
print(coverage_df[['pathway', 'total_genes', 'available_genes', 'coverage_percent']].head(10))

# Filter pathways with good coverage (>30% and >5 genes)
well_covered_pathways = coverage_df[
    (coverage_df['coverage_percent'] > 30) & 
    (coverage_df['available_genes'] > 5)
]['pathway'].tolist()

# Create filtered gene sets for analysis
filtered_metabolism_gene_sets = {
    pathway: [g for g in metabolism_gene_sets[pathway] if g in malignant_cells.var_names]
    for pathway in well_covered_pathways
}

# Extract all unique metabolism genes for NMF
all_metabolism_genes = set()
for pathway, genes in filtered_metabolism_gene_sets.items():
    all_metabolism_genes.update(genes)

all_metabolism_genes = list(all_metabolism_genes)
print(f"Total unique metabolism genes: {len(all_metabolism_genes)}")

# Check availability in malignant cells dataset
available_metabolism_genes = [gene for gene in all_metabolism_genes if gene in malignant_cells.var_names]
print(f"Available in dataset: {len(available_metabolism_genes)}")
print(f"Coverage: {len(available_metabolism_genes)/len(all_metabolism_genes)*100:.1f}%")
print(f"\n🎯 Ready for NMF with {len(available_metabolism_genes)} metabolism genes")

#%% Task 4: Robust Per-Patient NMF Analysis (Following CSCC Strategy)
print("\n🔬 Task 4: Robust Per-Patient Metabolism Program Discovery")
print("="*60)

## Use HVG genes for NMF    
if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

# adata = malignant_cells.copy()
# if adata.X.max() > 50:  # Check if data looks like raw counts
#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = 2000)
# available_metabolism_genes = adata.var_names[adata.var['highly_variable']].tolist()

# Run the robust NMF analysis
robust_results = run_robust_nmf_analysis(malignant_cells, available_metabolism_genes, result_path)

# Manually merge
updated_results, updated_cells = merge_meta_programs(
    malignant_cells=malignant_cells,
    robust_results=robust_results,
    source_mp=3,  # MP to merge (eliminate) # 7
    target_mp=2,  # MP to merge into (expand) # 1
    save_path=result_path,
    update_scoring=True
)

# Manually split
# updated_results, updated_cells = split_meta_program(
#     malignant_cells=malignant_cells,
#     robust_results=robust_results,
#     mp_to_split=4,
#     n_splits=2,
#     split_method='program_clustering',
#     save_path=result_path
# )

# Validate the merger
validation_stats = validate_merger(updated_results, updated_cells, save_path=result_path)

malignant_cells = updated_cells.copy()
robust_results = updated_results.copy()

# Load the robust results if not already loaded
if not 'robust_results' in locals():
    with open(os.path.join(result_path, "robust_metabolism_nmf_results.pkl"), 'rb') as f:
        robust_results = pickle.load(f)

# Save robust NMF results
print("\n💾 Exporting Results for R Statistical Analysis")
print("="*60)

export_dir = export_results_for_r_analysis(malignant_cells, robust_results, result_path)
print(f"\n🎉 Ready for R analysis!")
print(f"📁 Data exported to: {export_dir}")

# Save final results
malignant_cells.write(os.path.join(result_path, "malignant_cells_with_robust_mps.h5ad"))
# malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_cells_with_robust_mps.h5ad"))

print(f"\n🎉 Robust Meta-Program Analysis Completed!")
print(f"📁 Results saved to: {result_path}")

#%% Downstream analysis and visualizations
## Plot UMAP with Meta-Programs
print("\n📊 Plotting UMAP with Meta-Programs")
available_metabolism_genes = pd.read_csv(os.path.join(result_path, "NMF_results_for_R", "gene_names.csv"))["gene_name"].tolist()
create_mp_umap_visualization(malignant_cells, available_metabolism_genes, figurePath)

# Enhanced Clinical Analysis with Meta-Programs (move to R)
print("\n🔍 Enhanced Clinical Analysis with Meta-Programs")
print("="*70)
    
# Run enhanced analyses
print("\n💊 Meta-Program Distribution by Treatment Strategy")
treatment_mp_results = analyze_mp_distribution(malignant_cells, 'Treatment_Strategy', 
                                                'Meta-Programs by Treatment Strategy', figurePath)

print("\n🧬 Meta-Program Distribution by Microsatellite Status") 
msi_mp_results = analyze_mp_distribution(malignant_cells, 'Microsatellite_Status',
                                        'Meta-Programs by Microsatellite Status', figurePath)

adata = malignant_cells.copy()
task1_results = task1_intrinsic_resistance_msi_pretreatment(adata, result_path)
task2_results = task2_acquired_resistance_msi_nonpcr(adata, result_path)
task3_results = task3_msi_mss_mp_similarity(adata, result_path)

# Save all MP results
sample_mp_df = calculate_mp_fractions_per_sample(adata, mp_column='MP_assignment')
sample_mp_df.to_csv(f"{save_path}/all_sample_mp_fraction_data.csv", index=False)

print(f"\n✅ Complete analysis pipeline finished!")

# Other figures
# Figure 1: NMF Robustness Analysis
create_nmf_robustness_plots(robust_results, figurePath)

# Figure 2: MP Assignment Quality
create_mp_assignment_quality_plots(malignant_cells, figurePath)

# Figure 3: Hierarchical Clustering Dendrogram
create_clustering_dendrogram(robust_results, figurePath)

# Figure 4: Pathway Enrichment Analysis
create_pathway_enrichment_heatmap(robust_results, figurePath)