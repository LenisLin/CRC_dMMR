## Functions for pre-processing
import os
import gzip
import h5py

import scanpy as sc
import anndata as ad

import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix

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
            print(f"✓ Loaded from .h5 format: {h5_file}")
            
        elif study_name == 'GSE205506':
            # Format: 10X format (mtx + barcodes + features) - multiple samples
            adata = load_10x_multi_samples(data_path, files)
            print(f"✓ Loaded from 10X multi-sample format")
            
        elif study_name == 'GSE132257':
            # Format: Raw UMI count matrix + annotation
            count_file = [f for f in files if 'raw_UMI_count_matrix' in f][0]
            annotation_file = [f for f in files if 'cell_annotation' in f][0]
            adata = load_txt_format(data_path, count_file, annotation_file)
            print(f"✓ Loaded from txt format: {count_file}")
            
        elif study_name == 'GSE236581':
            # Format: Single 10X format (mtx + barcodes + features)
                """Load single 10X sample"""
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

                metadata_file = [f for f in files if 'metadata' in f][0]
                metadata = pd.read_csv(os.path.join(data_path, metadata_file), sep='\t')
                # Add metadata to adata.obs
                adata.obs = adata.obs.join(metadata.set_index(metadata.columns[0]), how='left')
                print(f"✓ Loaded from single 10X format with metadata: {metadata_file}")
            
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
            sample_id = file.split('_')[0]  # Extract sample ID
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
            # Read 10X data
            adata_sample = sc.read_10x_mtx(
                data_path,
                var_names='gene_symbols',
                cache=True,
                prefix=f"{sample_files['matrix'].split('_matrix')[0]}_"
            )
            adata_sample.var_names_make_unique()
            adata_sample.obs['sample_id'] = sample_id
            adatas.append(adata_sample)
            print(f"  ✓ Loaded sample {sample_id}: {adata_sample.shape}")
        except Exception as e:
            print(f"  ❌ Error loading sample {sample_id}: {str(e)}")
    
    # Concatenate all samples
    if adatas:
        adata_combined = ad.concat(adatas, axis=0, join='outer')
        adata_combined.obs_names_unique()
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
    
    # Create AnnData object
    adata = ad.AnnData(X=count_matrix.T.values)  # Transpose to cells x genes
    adata.obs_names = count_matrix.columns
    adata.var_names = count_matrix.index
    
    # Add annotation
    if annotation.shape[0] == adata.n_obs:
        # If annotation rows match cells
        annotation.index = adata.obs_names
        adata.obs = annotation
    else:
        print(f"⚠️  Annotation shape mismatch: {annotation.shape[0]} vs {adata.n_obs} cells")
        # Try to match by cell names
        if annotation.columns[0] in ['Cell', 'cell', 'barcode', 'Cell_name']:
            annotation.set_index(annotation.columns[0], inplace=True)
            adata.obs = adata.obs.join(annotation, how='left')
        else:
            print("Cannot match annotation to cells")
    
    return adata

def check_data_characteristics(adata, study_name):
    """Check and print data characteristics"""
    
    print(f"\n📊 Data Characteristics for {study_name}:")
    print(f"Shape: {adata.shape} (cells × genes)")
    
    # Check if data is raw counts
    print(f"\n🔍 Expression Data Check:")
    sample_values = adata.X[:100, :100] if hasattr(adata.X, 'toarray') else adata.X[:100, :100]
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

