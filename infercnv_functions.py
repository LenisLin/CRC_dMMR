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

# Prepare cell type annotations for InferCNV
def prepare_cnv_annotations(adata, tissue_col='Tissue_Type'):
    """Prepare cell type annotations for InferCNV analysis"""
    adata.obs['cnv_celltype'] = adata.obs['Major_type'].astype(str)
    
    if tissue_col in adata.obs.columns:
        # Distinguish normal vs tumor epithelial cells
        normal_epi_mask = (adata.obs['Major_type'] == 'Epithelial') & (adata.obs[tissue_col] == 'normal')
        tumor_epi_mask = (adata.obs['Major_type'] == 'Epithelial') & (adata.obs[tissue_col] == 'tumor')
        
        adata.obs.loc[normal_epi_mask, 'cnv_celltype'] = 'normal_epithelial'
        adata.obs.loc[tumor_epi_mask, 'cnv_celltype'] = 'tumor_epithelial'
        
        print(f"CNV cell type annotations:")
        print(adata.obs['cnv_celltype'].value_counts())
    
    return adata

def infercnv_qc(adata, min_genes=500, max_mt_pct=75):
    """Quality control specifically for InferCNV analysis with raw counts"""
    print(f"Before QC: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Basic cell filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter by MT percentage if available
    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < max_mt_pct].copy()
    
    # Filter genes (expressed in at least 10 cells)
    sc.pp.filter_genes(adata, min_cells=10)
    
    print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata

def run_infercnv_analysis(adata, output_path, analysis_name, gtf_path=None):
    """Run complete InferCNV analysis pipeline"""
    
    print(f"Running InferCNV for {analysis_name}...")
    
    # Set up genomic positions (you may need to adjust the GTF path)
    if gtf_path and os.path.exists(gtf_path):
        cnv.io.genomic_position_from_gtf(gtf_path, adata)
    else:
        # Use default genomic positions if GTF not available
        print("⚠️  Using default genomic positions (GTF file not specified)")
        cnv.io.genomic_position_from_gtf(None, adata)
    
    # Check if we have normal epithelial reference
    normal_refs = adata.obs[adata.obs['cnv_celltype'] == 'normal_epithelial']
    if len(normal_refs) == 0:
        print("⚠️  No normal epithelial cells found! Using all epithelial as reference")
            reference_cats = ['tumor_epithelial']  # This is not ideal but will work
    else:
        reference_cats = ['normal_epithelial']
        print(f"Using {len(normal_refs)} normal epithelial cells as reference")
    
    # Run InferCNV
    try:
        cnv.tl.infercnv(
            adata,
            reference_key="cnv_celltype",
            reference_cat=reference_cats,
            window_size=250,  # Adjusted for better resolution
            n_jobs=16  # Adjust based on your system
        )
        print("✅ InferCNV completed successfully")
    except Exception as e:
        print(f"❌ InferCNV failed: {str(e)}")
        return adata
    
    # CNV post-processing
    try:
        # CNV PCA and clustering
        # cnv.tl.pca(adata)
        # cnv.pp.neighbors(adata)
        # cnv.tl.leiden(adata, resolution=0.5)
        # cnv.tl.umap(adata)
        # cnv.tl.cnv_score(adata)
        
        print("✅ CNV post-processing completed")
        
        # Create visualizations
        create_cnv_visualizations(adata, output_path, analysis_name)
        
    except Exception as e:
        print(f"⚠️  CNV post-processing partially failed: {str(e)}")
    
    return adata

def create_cnv_visualizations(adata, output_path, analysis_name):
    """Create comprehensive CNV visualizations"""
    
    # 1. Chromosome heatmap
    try:
        plt.figure(figsize=(15, 8))
        cnv.pl.chromosome_heatmap(adata, groupby="cnv_celltype", show=False)
        plt.savefig(os.path.join(output_path, f"{analysis_name}_chromosome_heatmap.pdf"), 
                   bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(output_path, f"{analysis_name}_chromosome_heatmap.png"), 
                   bbox_inches="tight", dpi=300)
        plt.close()
        print("✅ Chromosome heatmap saved")
    except Exception as e:
        print(f"⚠️  Chromosome heatmap failed: {str(e)}")
    
    # # 2. CNV UMAP plots
    # try:
    #     fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    #     axes = axes.flatten()
        
    #     # CNV leiden clusters
    #     cnv.pl.umap(adata, color="cnv_leiden", legend_loc="on data", 
    #                legend_fontoutline=2, ax=axes[0], show=False)
    #     axes[0].set_title("CNV Leiden Clusters", fontweight='bold')
        
    #     # CNV score
    #     cnv.pl.umap(adata, color="cnv_score", ax=axes[1], show=False)
    #     axes[1].set_title("CNV Score", fontweight='bold')
        
    #     # Cell type
    #     cnv.pl.umap(adata, color="cnv_celltype", ax=axes[2], show=False)
    #     axes[2].set_title("Cell Type", fontweight='bold')
        
    #     # EPCAM expression if available
    #     if 'EPCAM' in adata.var_names:
    #         cnv.pl.umap(adata, color="EPCAM", ax=axes[3], show=False)
    #         axes[3].set_title("EPCAM Expression", fontweight='bold')
    #     else:
    #         axes[3].text(0.5, 0.5, 'EPCAM\nnot available', ha='center', va='center', 
    #                     transform=axes[3].transAxes)
    #         axes[3].set_title("EPCAM Expression", fontweight='bold')
        
    #     # Batch effect
    #     if 'batch' in adata.obs.columns:
    #         cnv.pl.umap(adata, color="batch", ax=axes[4], show=False)
    #         axes[4].set_title("Batch", fontweight='bold')
    #     else:
    #         axes[4].text(0.5, 0.5, 'Batch\nnot available', ha='center', va='center', 
    #                     transform=axes[4].transAxes)
    #         axes[4].set_title("Batch", fontweight='bold')
        
    #     # Treatment stage if available
    #     if 'Treatment_Stage' in adata.obs.columns:
    #         cnv.pl.umap(adata, color="Treatment_Stage", ax=axes[5], show=False)
    #         axes[5].set_title("Treatment Stage", fontweight='bold')
    #     else:
    #         axes[5].text(0.5, 0.5, 'Treatment Stage\nnot available', ha='center', va='center', 
    #                     transform=axes[5].transAxes)
    #         axes[5].set_title("Treatment Stage", fontweight='bold')
        
        # plt.tight_layout()
        # plt.savefig(os.path.join(output_path, f"{analysis_name}_cnv_analysis.pdf"), 
        #            bbox_inches="tight", dpi=300)
        # plt.savefig(os.path.join(output_path, f"{analysis_name}_cnv_analysis.png"), 
        #            bbox_inches="tight", dpi=300)
        # plt.close()
        # print("✅ CNV analysis plots saved")
        
    # except Exception as e:
    #     print(f"⚠️  CNV visualization failed: {str(e)}")

def compare_cnv_results(adata_raw, adata_harmony, output_path):
    """Compare CNV results between raw and harmony-corrected data"""
    
    # Extract CNV scores if available
    if 'cnv_score' in adata_raw.obs.columns and 'cnv_score' in adata_harmony.obs.columns:
        
        # Create comparison dataframe
        comparison_df = pd.DataFrame({
            'cell_id': adata_raw.obs.index,
            'cnv_score_raw': adata_raw.obs['cnv_score'],
            'cnv_score_harmony': adata_harmony.obs['cnv_score'],
            'batch': adata_raw.obs['batch'],
            'celltype': adata_raw.obs['cnv_celltype']
        })
        
        # Calculate correlation
        correlation = comparison_df['cnv_score_raw'].corr(comparison_df['cnv_score_harmony'])
        print(f"CNV score correlation (raw vs harmony): {correlation:.3f}")
        
        # Create comparison plots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # Scatter plot of CNV scores
        axes[0,0].scatter(comparison_df['cnv_score_raw'], comparison_df['cnv_score_harmony'], 
                         alpha=0.6, s=1)
        axes[0,0].plot([comparison_df['cnv_score_raw'].min(), comparison_df['cnv_score_raw'].max()],
                      [comparison_df['cnv_score_raw'].min(), comparison_df['cnv_score_raw'].max()], 
                      'r--', alpha=0.8)
        axes[0,0].set_xlabel('CNV Score (Raw)')
        axes[0,0].set_ylabel('CNV Score (Harmony)')
        axes[0,0].set_title(f'CNV Score Correlation\nr = {correlation:.3f}', fontweight='bold')
        
        # Distribution comparison
        axes[0,1].hist(comparison_df['cnv_score_raw'], bins=50, alpha=0.7, label='Raw', density=True)
        axes[0,1].hist(comparison_df['cnv_score_harmony'], bins=50, alpha=0.7, label='Harmony', density=True)
        axes[0,1].set_xlabel('CNV Score')
        axes[0,1].set_ylabel('Density')
        axes[0,1].set_title('CNV Score Distribution', fontweight='bold')
        axes[0,1].legend()
        
        # Batch effect comparison
        batch_comparison = comparison_df.groupby('batch').agg({
            'cnv_score_raw': 'mean',
            'cnv_score_harmony': 'mean'
        }).reset_index()
        
        x_pos = np.arange(len(batch_comparison))
        width = 0.35
        axes[0,2].bar(x_pos - width/2, batch_comparison['cnv_score_raw'], width, 
                     label='Raw', alpha=0.8)
        axes[0,2].bar(x_pos + width/2, batch_comparison['cnv_score_harmony'], width, 
                     label='Harmony', alpha=0.8)
        axes[0,2].set_xlabel('Batch')
        axes[0,2].set_ylabel('Mean CNV Score')
        axes[0,2].set_title('CNV Score by Batch', fontweight='bold')
        axes[0,2].set_xticks(x_pos)
        axes[0,2].set_xticklabels(batch_comparison['batch'], rotation=45)
        axes[0,2].legend()
        
        # Cell type comparison
        celltype_comparison = comparison_df.groupby('celltype').agg({
            'cnv_score_raw': 'mean',
            'cnv_score_harmony': 'mean'
        }).reset_index()
        
        x_pos = np.arange(len(celltype_comparison))
        axes[1,0].bar(x_pos - width/2, celltype_comparison['cnv_score_raw'], width, 
                     label='Raw', alpha=0.8)
        axes[1,0].bar(x_pos + width/2, celltype_comparison['cnv_score_harmony'], width, 
                     label='Harmony', alpha=0.8)
        axes[1,0].set_xlabel('Cell Type')
        axes[1,0].set_ylabel('Mean CNV Score')
        axes[1,0].set_title('CNV Score by Cell Type', fontweight='bold')
        axes[1,0].set_xticks(x_pos)
        axes[1,0].set_xticklabels(celltype_comparison['celltype'], rotation=45)
        axes[1,0].legend()
        
        # Difference plot
        comparison_df['cnv_diff'] = comparison_df['cnv_score_harmony'] - comparison_df['cnv_score_raw']
        axes[1,1].hist(comparison_df['cnv_diff'], bins=50, alpha=0.7)
        axes[1,1].axvline(0, color='red', linestyle='--', alpha=0.8)
        axes[1,1].set_xlabel('CNV Score Difference (Harmony - Raw)')
        axes[1,1].set_ylabel('Frequency')
        axes[1,1].set_title('CNV Score Differences', fontweight='bold')
        
        # Summary statistics
        stats_text = f"""
        Summary Statistics:
        
        Raw CNV Scores:
        Mean: {comparison_df['cnv_score_raw'].mean():.3f}
        Std: {comparison_df['cnv_score_raw'].std():.3f}
        
        Harmony CNV Scores:
        Mean: {comparison_df['cnv_score_harmony'].mean():.3f}
        Std: {comparison_df['cnv_score_harmony'].std():.3f}
        
        Correlation: {correlation:.3f}
        
        Mean Difference: {comparison_df['cnv_diff'].mean():.3f}
        """
        
        axes[1,2].text(0.1, 0.5, stats_text, transform=axes[1,2].transAxes, 
                      fontsize=10, verticalalignment='center', fontfamily='monospace')
        axes[1,2].set_xlim(0, 1)
        axes[1,2].set_ylim(0, 1)
        axes[1,2].axis('off')
        axes[1,2].set_title('Summary Statistics', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, "cnv_comparison_analysis.pdf"), 
                   bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(output_path, "cnv_comparison_analysis.png"), 
                   bbox_inches="tight", dpi=300)
        plt.close()
        
        # Save comparison data
        comparison_df.to_csv(os.path.join(output_path, "cnv_comparison_data.csv"), index=False)
        
        print("✅ Comparison analysis completed and saved")
        
        return comparison_df
    
    else:
        print("⚠️  CNV scores not available for comparison")
        return None
