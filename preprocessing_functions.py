## Functions for pre-processing
import os
import gzip
from turtle import st
import h5py
import warnings

import scanpy as sc
import scanpy.external as sce
import anndata as ad
import scrublet as scr

import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
from sklearn.metrics import silhouette_score

import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

#%% Data Loading functions
def load_scrna_dataset(study_name, data_path):
    """
    Load scRNA-seq dataset from various formats and check data characteristics
    
    Parameters:
    -----------
    study_name : str
        Name of the study (e.g., 'GSE108989')
    data_path : str
        Path to the study data directory
    
    Returns:
    --------
    adata : anndata.AnnData
        Loaded and basic processed AnnData object
    """
    
    print(f"\n{'='*60}")
    print(f"Processing Study: {study_name}")
    print(f"{'='*60}")
    
    files = os.listdir(data_path)
    print(f"Files found: {files}")
    
    # Initialize variables
    adata = None
    
    try:
        if study_name == 'GSE108989':
            # Format: .h5 file
            h5_file = [f for f in files if f.endswith('.h5')][0]
            adata = sc.read_10x_h5(os.path.join(data_path, h5_file))
            adata.var_names_make_unique()
            print(f"✓ Loaded from .h5 format: {h5_file}")
            
        elif study_name == 'GSE205506':
            # Format: 10X format (mtx + barcodes + features) - multiple samples
            adata = load_10x_multi_samples(data_path, files)

            print(f"✓ Loaded from 10X multi-sample format")
            
        elif study_name == 'GSE236581':
            # Format: Single 10X format (mtx + barcodes + features)
            matrix_file = [f for f in files if f.endswith('counts.mtx.gz')][0]
            barcodes_file = [f for f in files if f.endswith('barcodes.tsv.gz')][0]
            features_file = [f for f in files if f.endswith('features.tsv.gz')][0]
            
            # Read matrix
            matrix = mmread(os.path.join(data_path, matrix_file)).T.tocsr()
            
            # Read barcodes
            with gzip.open(os.path.join(data_path, barcodes_file), 'rt') as f:
                barcodes = [line.strip() for line in f]
            
            # Read features
            features = pd.read_csv(os.path.join(data_path, features_file), sep='\t', header=None)
            features.columns = ['gene_id', 'gene_symbol', 'feature_type']
            
            # Create AnnData object
            adata = ad.AnnData(X=matrix)
            adata.obs_names = barcodes
            adata.var_names = features['gene_symbol'].values
            adata.var['gene_id'] = features['gene_id'].values
            adata.var['feature_type'] = features['feature_type'].values
            adata.var_names_make_unique()

            # Load metadata if available
            metadata_files = [f for f in files if 'metadata' in f]
            if metadata_files:
                metadata_file = metadata_files[0]
                metadata = pd.read_csv(os.path.join(data_path, metadata_file), sep='\t')
                # Add metadata to adata.obs
                adata.obs = adata.obs.join(metadata.set_index(metadata.columns[0]), how='left')
                print(f"✓ Loaded from single 10X format with metadata: {metadata_file}")
            else:
                print(f"✓ Loaded from single 10X format (no metadata file found)")
            
        elif study_name == 'GSE144735':
            # Format: Raw UMI count matrix + annotation  
            count_file = [f for f in files if 'raw_UMI_count_matrix' in f][0]
            annotation_file = [f for f in files if 'annotation' in f][0]
            adata = load_txt_format(data_path, count_file, annotation_file)
            print(f"✓ Loaded from txt format: {count_file}")
            
        elif study_name == 'GSE178341':
            # Format: .h5 file
            h5_file = [f for f in files if f.endswith('.h5')][0]
            adata = sc.read_10x_h5(os.path.join(data_path, h5_file))
            adata.var_names_make_unique()
            adata.obs["Sample_ID"] = adata.obs_names.str.split('-').str[0]
            print(f"✓ Loaded from .h5 format: {h5_file}")
            
        elif study_name == 'GSE132465':
            # Format: Raw UMI count matrix + annotation
            count_file = [f for f in files if 'raw_UMI_count_matrix' in f][0]
            annotation_file = [f for f in files if 'cell_annotation' in f][0]
            adata = load_txt_format(data_path, count_file, annotation_file)
            print(f"✓ Loaded from txt format: {count_file}")
            
        else:
            print(f"❌ Unknown study format: {study_name}")
            return None
            
        # Add study information
        adata.obs['study'] = study_name
        
        # Check and report data characteristics
        check_data_characteristics(adata, study_name)
        
        return adata
        
    except Exception as e:
        print(f"❌ Error loading {study_name}: {str(e)}")
        return None

def load_10x_multi_samples(data_path, files):
    """Load multiple 10X samples and combine them"""
    
    # Group files by sample
    samples = {}
    for file in files:
        if 'matrix.mtx' in file:
            sample_id = file.split('_')[0]  # Extract sample ID (GSM number)
            if sample_id not in samples:
                samples[sample_id] = {}
            samples[sample_id]['matrix'] = file
        elif 'barcodes.tsv' in file:
            sample_id = file.split('_')[0]
            if sample_id not in samples:
                samples[sample_id] = {}
            samples[sample_id]['barcodes'] = file
        elif 'features.tsv' in file:
            sample_id = file.split('_')[0]
            if sample_id not in samples:
                samples[sample_id] = {}
            samples[sample_id]['features'] = file
    
    print(f"Found {len(samples)} samples: {list(samples.keys())}")
    
    # Load each sample and combine
    adatas = []
    for sample_id, sample_files in samples.items():
        try:
            # Check if all required files are present
            if not all(key in sample_files for key in ['matrix', 'barcodes', 'features']):
                print(f"  ⚠️ Sample {sample_id} missing required files, skipping")
                continue
                
            # Read individual files manually for better control
            matrix_path = os.path.join(data_path, sample_files['matrix'])
            barcodes_path = os.path.join(data_path, sample_files['barcodes'])
            features_path = os.path.join(data_path, sample_files['features'])
            
            # Read matrix
            matrix = mmread(matrix_path).T.tocsr()
            
            # Read barcodes
            with gzip.open(barcodes_path, 'rt') as f:
                barcodes = [line.strip() for line in f]
            
            # Read features
            features = pd.read_csv(features_path, sep='\t', header=None, compression='gzip')
            features.columns = ['gene_id', 'gene_symbol', 'feature_type']
            
            # Create AnnData object
            adata_sample = ad.AnnData(X=matrix)
            adata_sample.obs_names = barcodes
            adata_sample.var_names = features['gene_symbol'].values
            adata_sample.var['gene_id'] = features['gene_id'].values
            adata_sample.var['feature_type'] = features['feature_type'].values
            
            # Make gene names unique and add sample info
            adata_sample.var_names_make_unique()
            adata_sample.obs['sample_id'] = sample_id
            
            adatas.append(adata_sample)
            print(f"  ✓ Loaded sample {sample_id}: {adata_sample.shape}")
            
        except Exception as e:
            print(f"  ❌ Error loading sample {sample_id}: {str(e)}")
    
    # Concatenate all samples
    if adatas:
        adata_combined = ad.concat(adatas, axis=0, join='outer', merge='unique')
        adata_combined.obs["Sample_ID"] = adata_combined.obs["sample_id"]
        return adata_combined
    else:
        raise ValueError("No samples loaded successfully")

def load_txt_format(data_path, count_file, annotation_file):
    """Load data from txt format"""
    
    # Read count matrix
    if count_file.endswith('.gz'):
        count_matrix = pd.read_csv(os.path.join(data_path, count_file), 
                                 sep='\t', index_col=0, compression='gzip')
    else:
        count_matrix = pd.read_csv(os.path.join(data_path, count_file), 
                                 sep='\t', index_col=0)
    
    # Read annotation
    annotation = pd.read_csv(os.path.join(data_path, annotation_file), sep='\t')
    annotation = annotation.iloc[:,:4]
    
    # Create AnnData object
    adata = ad.AnnData(X=count_matrix.T.values)  # Transpose to cells x genes
    adata.obs_names = count_matrix.columns
    adata.var_names = count_matrix.index
    adata.var_names_make_unique()
    
    # Add annotation
    if annotation.shape[0] == adata.n_obs:
        # If annotation rows match cells
        annotation.index = adata.obs_names
        adata.obs = pd.concat([adata.obs, annotation], axis=1)
    else:
        print(f"⚠️  Annotation shape mismatch: {annotation.shape[0]} vs {adata.n_obs} cells")
        # Try to match by cell names
        if annotation.columns[0] in ['Cell', 'cell', 'barcode', 'Cell_name']:
            annotation_indexed = annotation.set_index(annotation.columns[0])
            adata.obs = adata.obs.join(annotation_indexed, how='left')
        else:
            print("Cannot match annotation to cells")
    
    return adata

def check_data_characteristics(adata, study_name):
    """Check and print data characteristics"""
    
    print(f"\n📊 Data Characteristics for {study_name}:")
    print(f"Shape: {adata.shape} (cells × genes)")
    
    # Check if data is raw counts
    print(f"\n🔍 Expression Data Check:")
    sample_values = adata.X[:100, :100]
    if hasattr(sample_values, 'toarray'):
        sample_values = sample_values.toarray()
    
    # Check for non-integer values
    has_decimals = np.any(sample_values % 1 != 0)
    max_val = np.max(sample_values)
    min_val = np.min(sample_values)
    
    print(f"  Max value: {max_val}")
    print(f"  Min value: {min_val}")
    print(f"  Has decimal values: {has_decimals}")
    print(f"  Data type: {sample_values.dtype}")
    
    if has_decimals or max_val < 1:
        print("  ⚠️  WARNING: Data appears to be normalized/log-transformed, not raw counts!")
    else:
        print("  ✓ Data appears to be raw counts")
    
    # Check gene names
    print(f"\n🧬 Gene Names Check:")
    print(f"  Total genes: {adata.n_vars}")
    print(f"  Gene name examples: {list(adata.var_names[:10])}")
    
    # Check for gene ID vs symbol
    if any(gene.startswith('ENSG') for gene in adata.var_names[:100]):
        print("  📝 Gene names appear to be Ensembl IDs")
    elif any(len(gene) > 10 and '|' in gene for gene in adata.var_names[:10]):
        print("  📝 Gene names appear to be composite (ID|Symbol)")
    else:
        print("  📝 Gene names appear to be gene symbols")
    
    # Check clinical/metadata information
    print(f"\n📋 Clinical/Metadata Information:")
    print(f"  Metadata columns: {list(adata.obs.columns)}")
    print(f"  Sample metadata preview:")
    print(adata.obs.head())
    
    # Check for key clinical variables
    key_vars = ['response', 'treatment', 'timepoint', 'patient', 'sample', 'MSI', 'stage']
    found_vars = []
    for var in key_vars:
        matches = [col for col in adata.obs.columns if var.lower() in col.lower()]
        if matches:
            found_vars.extend(matches)
    
    if found_vars:
        print(f"  🎯 Potential clinical variables found: {found_vars}")
    else:
        print("  ⚠️  No obvious clinical variables found")
    
    print(f"\n{'='*60}")

#%% Data preprocessing functions
def integrate_clinical_metadata(adata, study_name, clinical_file):
    """
    Integrate clinical metadata for studies that have external clinical data
    """
    
    if not os.path.exists(clinical_file):
        print(f"  ⚠️  Clinical data file not found: {clinical_file}")
        return adata
    
    try:
        clinical_data = pd.read_csv(clinical_file)
        print(f"  📋 Loading clinical data: {clinical_data.shape}")
        print(f"  Clinical columns: {list(clinical_data.columns)}")
        
        # Check for Sample_ID column
        if 'Sample_ID' not in clinical_data.columns:
            print(f"  ⚠️  'Sample_ID' column not found in clinical data")
            print(f"  Available columns: {list(clinical_data.columns)}")
            return adata
        
        # Extract sample IDs from cell barcodes or existing column
        if 'Sample_ID' in adata.obs.columns:
            # Use existing sample_id column
            sample_ids = adata.obs['Sample_ID']
        else:
            # Try to extract from obs_names (assuming format: SAMPLE_BARCODE)
            sample_ids = adata.obs_names.str.split('_').str[0]
            adata.obs['Sample_ID'] = sample_ids
        
        print(f"  🔍 Unique samples in data: {sample_ids.unique()}")
        print(f"  🔍 Unique samples in clinical: {clinical_data['Sample_ID'].unique()}")
        
        # Match samples
        clinical_data_indexed = clinical_data.set_index('Sample_ID')
        
        # Add clinical data to obs
        for col in clinical_data.columns:
            if col != 'Sample_ID':
                adata.obs[col] = adata.obs['Sample_ID'].map(clinical_data_indexed[col])
        
        # Report matching success
        if len(clinical_data.columns) > 1:
            matched_cells = adata.obs[clinical_data.columns[1]].notna().sum()
            print(f"  ✅ Successfully matched {matched_cells}/{adata.n_obs} cells with clinical data")
        
        return adata
        
    except Exception as e:
        print(f"  ❌ Error integrating clinical data: {str(e)}")
        return adata

def detect_doublets_scrublet(adata, study_name):
    """
    Detect doublets using Scrublet per study
    """
    print(f"  🔍 Running Scrublet doublet detection for {study_name}...")
    
    try:
        # Initialize Scrublet
        scrub = scr.Scrublet(adata.X)
        
        # Run doublet detection
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2, 
            min_cells=3, 
            min_gene_variability_pctl=85,
            n_prin_comps=30
        )
        
        # Add results to adata
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
        
        n_doublets = predicted_doublets.sum()
        doublet_rate = n_doublets / len(predicted_doublets) * 100
        print(f"  📊 Detected {n_doublets} doublets ({doublet_rate:.1f}% of cells)")

        # sc.pp.scrublet(adata, batch_key="sampleid")
        
        return adata
        
    except Exception as e:
        print(f"  ❌ Error in doublet detection: {str(e)}")
        # Add empty columns if failed
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata

def calculate_qc_metrics(adata):
    """
    Calculate QC metrics including mitochondrial and ribosomal gene content
    """
    print("  📏 Calculating QC metrics...")
    
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Identify ribosomal genes
    adata.var['ribo'] = (adata.var_names.str.startswith('RPS') | 
                        adata.var_names.str.startswith('RPL'))
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt",'ribo'], percent_top=None, log1p=False, inplace=True)
    
    print(f"    • Mitochondrial genes: {adata.var['mt'].sum()}")
    print(f"    • Ribosomal genes: {adata.var['ribo'].sum()}")
    print(f"    • Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}")
    print(f"    • Mean ribosomal %: {adata.obs['pct_counts_ribo'].mean():.2f}")
    
    return adata

def apply_quality_filters(adata, study_name):
    """
    Apply quality control filters based on original study specifications
    """
    print(f"  🔧 Applying study-specific quality filters for {study_name}...")
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    if study_name == 'GSE236581':
        # GSE236581 criteria:
        # (1) 600-25,000 UMI counts
        # (2) ≥600 detected genes  
        # (3) ≤5% MT (or ≤70% for CD45- cells, but we'll use 5% as standard)
        
        print(f"    📋 Applying GSE236581 filters (600-25K UMI, ≥600 genes, ≤5% MT)")
        
        # Filter by UMI counts
        adata = adata[(adata.obs['total_counts'] >= 600) & 
                     (adata.obs['total_counts'] <= 25000), :].copy()
        cells_after_umi = adata.n_obs
        print(f"    • UMI filter: {n_cells_before - cells_after_umi} cells removed (kept 600-25K UMI)")
        
        # Filter by gene count
        adata = adata[adata.obs['n_genes_by_counts'] >= 600, :].copy()
        cells_after_genes = adata.n_obs
        print(f"    • Gene filter: {cells_after_umi - cells_after_genes} cells removed (kept ≥600 genes)")
        
        # Filter by mitochondrial percentage
        adata = adata[adata.obs['pct_counts_mt'] <= 5, :].copy()
        cells_after_mt = adata.n_obs
        print(f"    • MT filter: {cells_after_genes - cells_after_mt} cells removed (kept ≤5% MT)")
        
    elif study_name == 'GSE205506':
        # GSE205506 criteria:
        # (1) 400-25,000 UMI counts
        # (2) 500-5,000 detected genes
        # (3) MT filtering done after cell type identification (skip for now)
        # Gene filtering: genes in <3 cells removed
        
        print(f"    📋 Applying GSE205506 filters (400-25K UMI, 500-5K genes, gene freq ≥3 cells)")
        
        # Filter genes first (genes in <3 cells)
        sc.pp.filter_genes(adata, min_cells=3)
        genes_after_filter = adata.n_vars
        print(f"    • Gene frequency filter: {n_genes_before - genes_after_filter} genes removed (kept genes in ≥3 cells)")
        
        # Filter by UMI counts
        adata = adata[(adata.obs['total_counts'] >= 400) & 
                     (adata.obs['total_counts'] <= 25000), :].copy()
        cells_after_umi = adata.n_obs
        print(f"    • UMI filter: {n_cells_before - cells_after_umi} cells removed (kept 400-25K UMI)")
        
        # Filter by gene count
        adata = adata[(adata.obs['n_genes_by_counts'] >= 500) & 
                     (adata.obs['n_genes_by_counts'] <= 5000), :].copy()
        cells_after_genes = adata.n_obs
        print(f"    • Gene filter: {cells_after_umi - cells_after_genes} cells removed (kept 500-5K genes)")
        
        print(f"    • Note: MT filtering will be applied after cell type identification per original protocol")
        
    elif study_name == 'GSE178341':
        # GSE178341 criteria:
        # (1) ≥200 genes
        # (2) ≥1,000 reads (using total_counts as proxy)
        # (3) ≥500 UMIs
        # (4) ≤50% mitochondrial UMIs
        # Additional outlier detection (simplified here)
        
        print(f"    📋 Applying GSE178341 filters (≥200 genes, ≥1K reads, ≥500 UMI, ≤50% MT)")
        
        # Filter by gene count
        adata = adata[adata.obs['n_genes_by_counts'] >= 200, :].copy()
        cells_after_genes = adata.n_obs
        print(f"    • Gene filter: {n_cells_before - cells_after_genes} cells removed (kept ≥200 genes)")
        
        # Filter by read count (using total_counts as proxy for reads)
        adata = adata[adata.obs['total_counts'] >= 1000, :].copy()
        cells_after_reads = adata.n_obs
        print(f"    • Read filter: {cells_after_genes - cells_after_reads} cells removed (kept ≥1K reads)")
        
        # Filter by UMI count (additional filter on top of total_counts)
        adata = adata[adata.obs['total_counts'] >= 500, :].copy()
        cells_after_umi = adata.n_obs
        print(f"    • UMI filter: {cells_after_reads - cells_after_umi} cells removed (kept ≥500 UMI)")
        
        # Filter by mitochondrial percentage
        adata = adata[adata.obs['pct_counts_mt'] <= 50, :].copy()
        cells_after_mt = adata.n_obs
        print(f"    • MT filter: {cells_after_umi - cells_after_mt} cells removed (kept ≤50% MT)")
        
    elif study_name in ['GSE132465', 'GSE144735']:
        # SMC and KUL3 datasets criteria:
        # (1) ≥1,000 UMI counts
        # (2) 200-6,000 genes
        # (3) ≤20% mitochondrial expression
        
        print(f"    📋 Applying {study_name} filters (≥1K UMI, 200-6K genes, ≤20% MT)")
        
        # Filter by UMI counts
        adata = adata[adata.obs['total_counts'] >= 1000, :].copy()
        cells_after_umi = adata.n_obs
        print(f"    • UMI filter: {n_cells_before - cells_after_umi} cells removed (kept ≥1K UMI)")
        
        # Filter by gene count
        adata = adata[(adata.obs['n_genes_by_counts'] >= 200) & 
                     (adata.obs['n_genes_by_counts'] <= 6000), :].copy()
        cells_after_genes = adata.n_obs
        print(f"    • Gene filter: {cells_after_umi - cells_after_genes} cells removed (kept 200-6K genes)")
        
        # Filter by mitochondrial percentage
        adata = adata[adata.obs['pct_counts_mt'] <= 20, :].copy()
        cells_after_mt = adata.n_obs
        print(f"    • MT filter: {cells_after_genes - cells_after_mt} cells removed (kept ≤20% MT)")
        
    elif study_name == 'GSE108989':
        # No specific filtering criteria mentioned, use conservative defaults
        print(f"    📋 Applying default filters for {study_name} (≥200 genes, ≤20% MT)")
        
        # Filter by gene count (conservative)
        adata = adata[adata.obs['n_genes_by_counts'] >= 200, :].copy()
        cells_after_genes = adata.n_obs
        print(f"    • Gene filter: {n_cells_before - cells_after_genes} cells removed (kept ≥200 genes)")
        
        # Filter by mitochondrial percentage
        adata = adata[adata.obs['pct_counts_mt'] <= 20, :].copy()
        cells_after_mt = adata.n_obs
        print(f"    • MT filter: {cells_after_genes - cells_after_mt} cells removed (kept ≤20% MT)")
        
    else:
        # Fallback for unknown studies
        print(f"    ⚠️  Unknown study {study_name}, applying default filters")
        sc.pp.filter_cells(adata, min_genes=200)
        adata = adata[adata.obs['pct_counts_mt'] <= 20, :].copy()
    
    # Common gene filtering (remove MT and ribosomal genes for all studies)
    genes_before_filter = adata.n_vars
    adata = adata[:, ~adata.var['mt']].copy()
    adata = adata[:, ~adata.var['ribo']].copy()
    print(f"    • Removed {genes_before_filter - adata.n_vars} MT/ribosomal genes (common filter)")
    
    print(f"    ✅ Final shape after {study_name}-specific filtering: {adata.shape}")
    
    return adata

def preprocess_single_dataset(adata, study_name, clinical_data_path=None):
    """
    Complete preprocessing pipeline for a single dataset
    """
    print(f"\n🔄 Preprocessing {study_name}...")
    
    # Step 1: Integrate clinical metadata
    if clinical_data_path:
        adata = integrate_clinical_metadata(adata, study_name, clinical_file = clinical_data_path)
    
    # Step 2: Calculate QC metrics
    adata = calculate_qc_metrics(adata)

    # Step 3: Apply quality filters
    adata = apply_quality_filters(adata, study_name)
        
    # Step 4: Detect doublets
    adata = detect_doublets_scrublet(adata, study_name)

    # Step 5: Basic preprocessing for downstream analysis
    print("  🔬 Running basic preprocessing...")
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    # Store raw data
    adata.raw = adata
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    print(f"    • Found {adata.var['highly_variable'].sum()} highly variable genes")
    
    # Add preprocessing info
    adata.uns['preprocessing_info'] = {
        'study': study_name,
        'min_genes_filter': 300,
        'max_mt_pct_filter': 20,
        'doublets_detected': adata.obs['predicted_doublet'].sum(),
        'n_hvgs': adata.var['highly_variable'].sum()
    }
    
    print(f"  ✅ Preprocessing completed for {study_name}")
    return adata

#%% Preparation for integration
def merge_datasets_for_integration(preprocessed_datasets):
    """
    Merge all preprocessed datasets into a single AnnData object for integration
    """
    print(f"\n{'='*60}")
    print("🔗 Merging datasets for integration")
    print(f"{'='*60}")
    
    # Collect all datasets
    adata_list = []
    for study_name, adata in preprocessed_datasets.items():
        print(f"  Adding {study_name}: {adata.shape}")
        adata_list.append(adata)
    
    # Merge datasets
    print("  🔄 Concatenating datasets...")
    merged_adata = ad.concat(adata_list, axis=0, join='outer', merge='unique')
    
    # Make observation names unique
    merged_adata.obs_names_make_unique()
    
    # Add batch information for integration
    merged_adata.obs['batch'] = merged_adata.obs["study"].astype('category')
    
    print(f"  ✅ Merged dataset shape: {merged_adata.shape}")
    print(f"  📊 Studies included: {merged_adata.obs['batch'].value_counts()}")
    
    return merged_adata


#%% Batch correction methods
def downsample_stratified(adata, n_sample=50000, batch_key='batch', 
                         quality_weight=True, random_state=42):
    """
    Perform stratified downsampling with optional quality weighting
    
    Parameters:
    -----------
    adata : AnnData
        Input data
    n_sample : int
        Target number of cells to sample
    batch_key : str
        Column name for batch stratification
    quality_weight : bool
        Whether to weight sampling by cell quality
    """
    print(f"🎯 Downsampling from {adata.n_obs:,} to {n_sample:,} cells...")
    
    np.random.seed(random_state)
    
    # Get batch information
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)
    
    # Calculate cells per batch (proportional to original)
    batch_counts = adata.obs[batch_key].value_counts()
    batch_proportions = batch_counts / batch_counts.sum()
    
    # Ensure minimum representation per batch
    min_per_batch = max(500, n_sample // (n_batches * 4))  # At least 500 or 1/4 of average
    
    sampled_indices = []
    
    for batch in batches:
        batch_mask = adata.obs[batch_key] == batch
        batch_cells = np.where(batch_mask)[0]
        
        # Calculate target for this batch
        target_batch = max(min_per_batch, int(n_sample * batch_proportions[batch]))
        target_batch = min(target_batch, len(batch_cells))  # Don't exceed available
        
        if quality_weight and 'n_genes_by_counts' in adata.obs.columns:
            # Weight by cell quality (higher gene count = higher probability)
            quality_scores = adata.obs.loc[batch_mask, 'n_genes_by_counts'].values
            # Normalize to probabilities
            weights = quality_scores / quality_scores.sum()
            
            # Sample with replacement=False
            batch_sample = np.random.choice(
                batch_cells, 
                size=target_batch, 
                replace=False, 
                p=weights
            )
        else:
            # Random sampling
            batch_sample = np.random.choice(
                batch_cells, 
                size=target_batch, 
                replace=False
            )
        
        sampled_indices.extend(batch_sample)
        print(f"  • {batch}: {len(batch_sample):,} cells (original: {len(batch_cells):,})")
    
    # Create downsampled data
    adata_downsampled = adata[sampled_indices].copy()
    
    print(f"  ✅ Downsampled to {adata_downsampled.n_obs:,} cells")
    print(f"  📊 Batch distribution: {adata_downsampled.obs[batch_key].value_counts().to_dict()}")
    
    return adata_downsampled

def prepare_for_integration(adata, n_top_genes=2000, downsample_n=None, 
                          batch_key='batch', quality_weight=True):
    """
    Prepare data for batch correction with optional downsampling
    """
    print("🔬 Preparing data for integration...")
    print(f"  Initial shape: {adata.shape}")
    
    # Downsample if requested
    if downsample_n is not None and adata.n_obs > downsample_n:
        adata = downsample_stratified(
            adata, 
            n_sample=downsample_n, 
            batch_key=batch_key,
            quality_weight=quality_weight
        )
    else:
        print(f"  No downsampling needed (target: {downsample_n}, actual: {adata.n_obs})")
    
    # Store raw data if not already done
    if adata.raw is None:
        adata.raw = adata
    
    # Find highly variable genes across all datasets
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=n_top_genes, 
        batch_key=batch_key,
        subset=True
    )
    
    print(f"  • Selected {adata.n_vars} highly variable genes")
    
    # Scale data for integration
    sc.pp.scale(adata, max_value=10)
    
    # Compute PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    
    print(f"  ✅ Data prepared: {adata.shape}")
    return adata

def apply_bbknn_integration(adata, batch_key='batch'):
    """
    Apply BBKNN batch correction
    """
    print(f"\n🔗 Applying BBKNN integration on {adata.n_obs:,} cells...")
    
    adata_bbknn = adata.copy()
    
    try:
        # Apply BBKNN
        sc.external.pp.bbknn(
            adata_bbknn, 
            batch_key=batch_key,
            neighbors_within_batch=3,
            n_pcs=50
        )
        
        # Compute UMAP
        sc.tl.umap(adata_bbknn)
        
        print("  ✅ BBKNN integration completed")
        return adata_bbknn
        
    except Exception as e:
        print(f"  ❌ BBKNN integration failed: {str(e)}")
        return None

def apply_harmony_integration(adata, batch_key='batch'):
    """
    Apply Harmony batch correction
    """
    print(f"\n🎵 Applying Harmony integration on {adata.n_obs:,} cells...")
    
    adata_harmony = adata.copy()
    
    try:
        # Apply Harmony
        sc.external.pp.harmony_integrate(
            adata_harmony, 
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_pca_harmony'
        )
        
        # Compute UMAP on Harmony-corrected data
        sc.pp.neighbors(adata_harmony, use_rep='X_pca_harmony', n_neighbors=15)
        sc.tl.umap(adata_harmony)
        
        print("  ✅ Harmony integration completed")
        return adata_harmony
        
    except Exception as e:
        print(f"  ❌ Harmony integration failed: {str(e)}")
        return None

def apply_mnn_integration(adata, batch_key='batch'):
    """
    Apply MNN (Mutual Nearest Neighbors) batch correction
    """
    print(f"\n🔗 Applying MNN integration on {adata.n_obs:,} cells...")
    
    adata_mnn = adata.copy()
    
    try:
        # Apply MNN correction
        sc.external.pp.mnn_correct(
            adata_mnn,
            batch_key=batch_key,
            save_raw=True,  # Save raw data
            n_jobs=1,  # Use single core to avoid memory issues
            batch_categories=None  # Auto-detect batch categories
        )
        
        # The corrected data is stored in adata.X
        # Compute PCA on corrected data
        sc.tl.pca(adata_mnn, svd_solver='arpack', n_comps=50)
        
        # Compute neighbors and UMAP
        sc.pp.neighbors(adata_mnn, n_neighbors=15, n_pcs=50)
        sc.tl.umap(adata_mnn)
        
        print("  ✅ MNN integration completed")
        return adata_mnn
        
    except Exception as e:
        print(f"  ❌ MNN integration failed: {str(e)}")
        return None

def apply_scanorama_integration(adata, batch_key='batch'):
    """
    Apply Scanorama batch correction
    """
    print(f"\n🔄 Applying Scanorama integration on {adata.n_obs:,} cells...")
    
    adata_scanorama = adata.copy()
    
    try:        
        # Apply Scanorama
        sc.external.pp.scanorama_integrate(
            adata_scanorama, 
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_scanorama'
        )  
        
        # Compute neighbors and UMAP on corrected data
        sc.pp.neighbors(adata_scanorama, use_rep='X_scanorama', n_neighbors=15)
        sc.tl.umap(adata_scanorama)
        
        print("  ✅ Scanorama integration completed")
        return adata_scanorama
        
    except Exception as e:
        print(f"  ❌ Scanorama integration failed: {str(e)}")
        return None

def compute_integration_metrics(adata, batch_key='batch', n_sample=5000):
    """
    Compute metrics to evaluate integration quality
    """
    # Sample data if too large for metrics computation
    adata_metrics = adata.copy()
    if adata_metrics.n_obs > n_sample:
        print(f"    📊 Subsampling to {n_sample} cells for metrics computation")
        sc.pp.subsample(adata_metrics, n_obs=n_sample)
    
    # Get batch labels
    batch_labels = adata_metrics.obs[batch_key].astype('category').cat.codes
    
    # Compute silhouette score (lower is better for batch correction)
    try:
        sil_score = silhouette_score(adata_metrics.obsm['X_umap'], batch_labels)
        return {'silhouette_batch': sil_score}
    except:
        return {'silhouette_batch': np.nan}

def visualize_integration_results(adata_dict, batch_key='batch', figsize=(20, 12), 
                                n_cells_plot=10000):
    """
    Create comprehensive visualization of integration results
    """
    print(f"\n📊 Creating integration comparison plots (max {n_cells_plot:,} cells per plot)...")
    
    n_methods = len(adata_dict)
    fig, axes = plt.subplots(2, n_methods, figsize=figsize)
    
    if n_methods == 1:
        axes = axes.reshape(2, 1)
    
    # Plot batch effects
    for i, (method, adata) in enumerate(adata_dict.items()):
        if adata is not None:
            # Sample cells for visualization if too many
            if adata.n_obs > n_cells_plot:
                print(f"    📊 Subsampling {method} to {n_cells_plot:,} cells for visualization")
                adata_vis = adata[np.random.choice(adata.n_obs, n_cells_plot, replace=False)]
            else:
                adata_vis = adata
            
            # Batch plot
            sc.pl.umap(
                adata_vis, 
                color=batch_key, 
                ax=axes[0, i], 
                show=False, 
                frameon=False,
                title=f'{method} - Batch ({adata_vis.n_obs:,} cells)'
            )
                
            # Study plot  
            sc.pl.umap(
                adata_vis, 
                color='study', 
                ax=axes[1, i], 
                show=False, 
                frameon=False,
                title=f'{method} - Study ({adata_vis.n_obs:,} cells)'
            )
        else:
            axes[0, i].text(0.5, 0.5, f'{method}\nFailed', 
                          ha='center', va='center', transform=axes[0, i].transAxes)
            axes[1, i].text(0.5, 0.5, f'{method}\nFailed', 
                          ha='center', va='center', transform=axes[1, i].transAxes)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def compare_all_integration_methods(adata_prepared, batch_key='batch'):
    """
    Apply all integration methods and compare results
    """
    print(f"\n{'='*80}")
    print(f"🚀 Comprehensive Batch Correction Comparison ({adata_prepared.n_obs:,} cells)")
    print(f"{'='*80}")
    
    # Store original (uncorrected) data with UMAP
    print(f"\n📍 Computing uncorrected UMAP on {adata_prepared.n_obs:,} cells...")
    adata_uncorrected = adata_prepared.copy()
    sc.pp.neighbors(adata_uncorrected, n_neighbors=15)
    sc.tl.umap(adata_uncorrected)
    
    # Dictionary to store results
    results = {'Uncorrected': adata_uncorrected}
    
    # Apply each method
    methods = [
        ('Harmony', apply_harmony_integration),
        ('BBKNN', apply_bbknn_integration),
        ('MNN', apply_mnn_integration),
        ('Scanorama', apply_scanorama_integration)
    ]
    
    for method_name, method_func in methods:
        try:
            result = method_func(adata_prepared.copy(), batch_key=batch_key)
            results[method_name] = result
        except Exception as e:
            print(f"❌ {method_name} failed: {str(e)}")
            results[method_name] = None
    
    # Compute metrics
    print(f"\n📊 Computing integration metrics...")
    metrics_df = pd.DataFrame()
    
    for method_name, adata_result in results.items():
        if adata_result is not None:
            metrics = compute_integration_metrics(adata_result.copy(), batch_key=batch_key)
            metrics['method'] = method_name
            metrics['n_cells'] = adata_result.n_obs
            metrics['n_genes'] = adata_result.n_vars
            metrics_df = pd.concat([metrics_df, pd.DataFrame([metrics])], ignore_index=True)
    
    # Print metrics
    print("\n📈 Integration Quality Metrics:")
    print("=" * 60)
    print(metrics_df[['method', 'silhouette_batch', 'n_cells', 'n_genes']].round(3))
    print("\nNote: Lower silhouette_batch score indicates better batch mixing")
    
    # Create visualizations
    fig = visualize_integration_results(results, batch_key=batch_key)
    
    # Print summary
    print(f"\n{'='*80}")
    print("🎯 Integration Summary:")
    successful_methods = [k for k, v in results.items() if v is not None]
    failed_methods = [k for k, v in results.items() if v is None]
    
    print(f"✅ Successful methods: {successful_methods}")
    if failed_methods:
        print(f"❌ Failed methods: {failed_methods}")
    
    print(f"{'='*80}")
    
    return results, metrics_df, fig

def save_integration_results(results, save_path="./integration_results/"):
    """
    Save integration results with downsampling info
    """
    import os
    os.makedirs(save_path, exist_ok=True)
    
    for method_name, adata in results.items():
        if adata is not None:
            filename = f"{save_path}/{method_name.lower()}_integrated_n{adata.n_obs}.h5ad"
            adata.write(filename)
            print(f"💾 Saved {method_name} results ({adata.n_obs:,} cells) to {filename}")

# Two-phase workflow functions
def quick_method_comparison(merged_adata, batch_key='batch', n_sample=20000, n_top_genes=2000):
    """
    Phase 1: Quick method comparison with small sample
    """
    print(f"\n🚀 PHASE 1: Quick Method Comparison (n={n_sample:,})")
    print("="*60)
    
    # Prepare data with downsampling
    adata_small = prepare_for_integration(
        merged_adata, 
        n_top_genes=n_top_genes,
        downsample_n=n_sample,
        batch_key=batch_key
    )
    
    # Compare methods
    results, metrics_df, fig = compare_all_integration_methods(adata_small, batch_key=batch_key)
    
    return results, metrics_df, fig