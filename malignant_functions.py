# Comprehensive Clinical Analysis of NMF Clusters
import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform, pdist

from statsmodels.stats.multitest import multipletests

import warnings
warnings.filterwarnings('ignore')

import gc
import pickle
import json
from collections import Counter, defaultdict

import gseapy as gp
from gseapy import Msigdb
from sklearn.decomposition import NMF
from itertools import combinations
import scipy.io
import scipy.sparse

def get_metabolism_pathways():
    """Enhanced metabolism pathway collection"""
    
    metabolism_gene_sets = {}
    
    # Enhanced metabolism keywords
    metabolism_keywords = [
        'metabolism', 'metabolic', 'glycolysis', 'glycolytic',
        'oxidative phosphorylation', 'electron transport', 'atp synthesis',
        'fatty acid', 'lipid', 'cholesterol', 'steroid',
        'amino acid', 'protein', 'nucleotide', 'purine', 'pyrimidine',
        'tca cycle', 'citrate cycle', 'krebs cycle',
        'pentose phosphate', 'hexosamine', 'fructose',
        'gluconeogenesis', 'glucose', 'insulin signaling',
        'autophagy', 'lysosome', 'peroxisome',
        'mitochondrial', 'ribosome', 'translation'
    ]
    
    # 1. KEGG Enhanced
    # try:
    #     kegg = gp.get_library(name='KEGG_2021_Human', organism='Human')
    #     kegg_metabolism = {k: v for k, v in kegg.items() 
    #                       if any(keyword in k.lower() for keyword in metabolism_keywords)}
    #     metabolism_gene_sets.update({f"KEGG_{k}": v for k, v in kegg_metabolism.items()})
    # except:
    #     pass
    
    # # 2. Reactome Enhanced
    # try:
    #     reactome = gp.get_library(name='Reactome_2022', organism='Human')
    #     reactome_metabolism = {k: v for k, v in reactome.items() 
    #                           if any(keyword in k.lower() for keyword in metabolism_keywords)}
    #     metabolism_gene_sets.update({f"REACTOME_{k}": v for k, v in reactome_metabolism.items()})
    # except:
    #     pass

    # 3. Hallmark Metabolism Pathways
    try:
        msig = Msigdb() # hallmark gene sets
        # msig.list_dbver() # list msigdb version you wanna query
        # msig.list_category(dbver="2025.1.Hs") # list categories given dbver.
        category = 'h.all'  # Hallmark gene sets
        dbver = "2025.1.Hs"  # Specify the version
        print(f"Fetching Hallmark gene sets from {category} with version {dbver}...")

        hallmark = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
        hallmark_metabolism = {k: v for k, v in hallmark.items() 
                              if any(keyword in k.lower() for keyword in metabolism_keywords)}
        
        metabolism_gene_sets.update({f"HALLMARK_{k}": v for k, v in hallmark_metabolism.items()})
        print(f"✅ Hallmark metabolism pathways: {len(hallmark_metabolism)}")
        
    except Exception as e:
        print(f"⚠️ Hallmark access failed: {e}")

    
    # # 4. GO Biological Process - Metabolism
    # try:
    #     go_bp = gp.get_library(name='GO_Biological_Process_2023', organism='Human')
    #     go_metabolism = {k: v for k, v in go_bp.items() 
    #                     if any(keyword in k.lower() for keyword in metabolism_keywords)}
    #     # Limit GO terms to prevent overload
    #     go_metabolism_top = dict(list(go_metabolism.items())[:100])
    #     metabolism_gene_sets.update({f"GO_BP_{k}": v for k, v in go_metabolism_top.items()})
    # except:
    #     pass
    
    return metabolism_gene_sets

def get_oncogenic_pathways():
    """Get oncogenic and cancer hallmark pathways"""
    
    oncogenic_gene_sets = {}
    
    # Cancer hallmark keywords
    cancer_keywords = [
        'cancer', 'tumor', 'oncogene', 'tumor suppressor',
        'cell cycle', 'apoptosis', 'proliferation', 'growth',
        'dna repair', 'dna damage', 'genome instability',
        'angiogenesis', 'invasion', 'metastasis',
        'emt', 'epithelial mesenchymal', 'stemness',
        'p53', 'rb', 'myc', 'ras', 'pi3k', 'akt', 'mtor',
        'wnt', 'notch', 'hedgehog', 'tgf', 'nf-kb'
    ]
    
    # 1. Hallmark Gene Sets
    try:
        hallmark = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
        oncogenic_gene_sets.update({f"HALLMARK_{k}": v for k, v in hallmark.items()})
    except:
        pass
    
    # # 2. Oncogenic Signatures
    # try:
    #     oncogenic = gp.get_library(name='MSigDB_Oncogenic_Signatures', organism='Human')
    #     oncogenic_gene_sets.update({f"ONCOGENIC_{k}": v for k, v in oncogenic.items()})
    # except:
    #     pass
    
    # # 3. KEGG Cancer Pathways
    # try:
    #     kegg = gp.get_library(name='KEGG_2021_Human', organism='Human')
    #     kegg_cancer = {k: v for k, v in kegg.items() 
    #                   if any(keyword in k.lower() for keyword in cancer_keywords)}
    #     oncogenic_gene_sets.update({f"KEGG_CANCER_{k}": v for k, v in kegg_cancer.items()})
    # except:
    #     pass
    
    return oncogenic_gene_sets

def get_immune_microenvironment_pathways():
    """Get immune and microenvironment pathways"""
    
    immune_gene_sets = {}
    
    immune_keywords = [
        'immune', 'immunity', 'immunological',
        'interferon', 'interleukin', 'cytokine', 'chemokine',
        'inflammatory', 'inflammation', 'response',
        't cell', 'b cell', 'nk cell', 'macrophage', 'dendritic',
        'antigen presentation', 'mhc', 'hla',
        'complement', 'toll like', 'innate immunity',
        'adaptive immunity', 'antibody', 'immunoglobulin'
    ]
    
    try:
        # Reactome immune pathways
        reactome = gp.get_library(name='Reactome_2022', organism='Human')
        reactome_immune = {k: v for k, v in reactome.items() 
                          if any(keyword in k.lower() for keyword in immune_keywords)}
        immune_gene_sets.update({f"REACTOME_IMMUNE_{k}": v for k, v in reactome_immune.items()})
        
        # GO immune process
        go_bp = gp.get_library(name='GO_Biological_Process_2023', organism='Human')
        go_immune = {k: v for k, v in go_bp.items() 
                    if any(keyword in k.lower() for keyword in immune_keywords)}
        go_immune_top = dict(list(go_immune.items())[:50])
        immune_gene_sets.update({f"GO_IMMUNE_{k}": v for k, v in go_immune_top.items()})
        
    except:
        pass
    
    return immune_gene_sets

def get_stress_resistance_pathways():
    """Get stress response and treatment resistance pathways"""
    
    stress_gene_sets = {}
    
    stress_keywords = [
        'stress', 'response', 'resistance', 'drug resistance',
        'hypoxia', 'oxidative stress', 'er stress', 'unfolded protein',
        'heat shock', 'chaperone', 'proteostasis',
        'autophagy', 'apoptosis', 'necroptosis', 'ferroptosis',
        'dna damage response', 'checkpoint', 'repair',
        'chemotherapy', 'radiation', 'therapy resistance'
    ]
    
    try:
        # Hallmark stress pathways
        hallmark = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
        hallmark_stress = {k: v for k, v in hallmark.items() 
                          if any(keyword in k.lower() for keyword in stress_keywords)}
        stress_gene_sets.update({f"HALLMARK_STRESS_{k}": v for k, v in hallmark_stress.items()})
        
        # # Reactome stress pathways
        # reactome = gp.get_library(name='Reactome_2022', organism='Human')
        # reactome_stress = {k: v for k, v in reactome.items() 
        #                   if any(keyword in k.lower() for keyword in stress_keywords)}
        # stress_gene_sets.update({f"REACTOME_STRESS_{k}": v for k, v in reactome_stress.items()})
        
    except:
        pass
    
    return stress_gene_sets

def get_highly_variable_genes(adata, candidate_genes, n_top_genes=3000):
    """Get highly variable genes from candidates"""
    
    if 'highly_variable' not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    
    hvg_genes = adata.var_names[adata.var.highly_variable].tolist()
    hvg_candidates = [gene for gene in hvg_genes if gene in candidate_genes]
    
    return hvg_candidates

def get_multi_pathway_genes(gene_sets, candidate_genes, min_pathways=2):
    """Get genes that appear in multiple pathways"""
    
    gene_pathway_count = defaultdict(int)
    
    for pathway, genes in gene_sets.items():
        for gene in genes:
            if gene in candidate_genes:
                gene_pathway_count[gene] += 1
    
    multi_pathway_genes = [gene for gene, count in gene_pathway_count.items() 
                          if count >= min_pathways]
    
    return multi_pathway_genes

def get_cancer_associated_genes(candidate_genes):
    """Get known cancer-associated genes"""
    
    # Common cancer genes (you can expand this list)
    cancer_genes = {
        'TP53', 'KRAS', 'PIK3CA', 'APC', 'PTEN', 'EGFR', 'MYC', 'RB1',
        'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'MLH1', 'MSH2', 'MSH6', 'PMS2',
        'VHL', 'NF1', 'NF2', 'CDKN2A', 'CDK4', 'MDM2', 'ERBB2', 'MET',
        'ALK', 'ROS1', 'RET', 'NTRK1', 'BRAF', 'NRAS', 'HRAS', 'KIT',
        'PDGFRA', 'FLT3', 'IDH1', 'IDH2', 'TERT', 'ARID1A', 'CTNNB1'
    }
    
    cancer_candidates = [gene for gene in candidate_genes if gene in cancer_genes]
    
    return cancer_candidates

def calculate_gene_overlap(genes1, genes2):
    """Calculate Jaccard similarity between two gene sets"""
    set1, set2 = set(genes1), set(genes2)
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union > 0 else 0

def get_overlap_count(genes1, genes2):
    """Get number of overlapping genes"""
    return len(set(genes1) & set(genes2))

def run_robust_nmf_analysis(malignant_cells, available_metabolism_genes, save_path):
    """
    Implement robust NMF strategy following CSCC paper methodology
    """
    
    print("🎯 Starting robust NMF analysis...")
    
    # Parameters following the paper
    K_RANGE = range(3, 12)  # Try: range(3, 8) or range(4, 10)
    TOP_GENES_PER_PROGRAM = 35  # Top 25 genes per program
    MIN_CELLS_PER_PATIENT = 150  # Minimum 150 cells per patient
    MAX_CLUSTER_NUM = 7
    
    # Robustness criteria
    SAME_PATIENT_OVERLAP_THRESHOLD = 0.7  # 70% overlap (35/50 genes)
    CROSS_PATIENT_OVERLAP_THRESHOLD = 0.3  # 30% overlap (10/50 genes)
    MAX_INTRA_PATIENT_OVERLAP = 0.3  # <20% overlap within same patient
    
    # Group cells by patient/study
    patient_groups = malignant_cells.obs['Sample_ID'].unique()
    print(f"Found {len(patient_groups)} unique samples")
    
    # Filter patients with sufficient cells
    valid_patients = []
    for patient in patient_groups:
        patient_cells = malignant_cells[malignant_cells.obs['Sample_ID'] == patient]
        if patient_cells.n_obs >= MIN_CELLS_PER_PATIENT:
            valid_patients.append(patient)
    
    print(f"Samples with ≥{MIN_CELLS_PER_PATIENT} cells: {len(valid_patients)}")
    
    # Step 1: Run NMF for each patient and each K
    print("\n📊 Step 1: Running NMF for each patient...")
    all_programs = []
    patient_programs = {}
    
    for patient_id in valid_patients:
        print(f"  Processing sample {patient_id}...")
        
        # Subset patient cells
        patient_cells = malignant_cells[malignant_cells.obs['Sample_ID'] == patient_id]
        patient_metabolism = patient_cells[:, available_metabolism_genes].copy()
        
        # Preprocessing
        sc.pp.normalize_total(patient_metabolism, target_sum=1e4)
        if scipy.sparse.issparse(patient_metabolism.X):
            X_patient = patient_metabolism.X.toarray()
        else:
            X_patient = patient_metabolism.X.copy()
        
        X_patient = np.maximum(X_patient, 0)  # Ensure non-negative
        
        patient_programs[patient_id] = {}
        
        # Run NMF for different K values
        for k in K_RANGE:
            print(f"    K={k}: {patient_metabolism.n_obs} cells, {patient_metabolism.n_vars} genes")
            
            # Run NMF
            nmf_model = NMF(n_components=k, random_state=619, max_iter=2000)
            W = nmf_model.fit_transform(X_patient)
            H = nmf_model.components_
            
            # Extract top genes for each program
            k_programs = []
            for factor_idx in range(k):
                # Get top 50 genes for this factor
                factor_scores = H[factor_idx, :]
                top_gene_indices = np.argsort(factor_scores)[-TOP_GENES_PER_PROGRAM:][::-1]
                top_genes = [available_metabolism_genes[i] for i in top_gene_indices]
                
                program = {
                    'patient_id': patient_id,
                    'k_value': k,
                    'factor_idx': factor_idx,
                    'genes': top_genes,
                    'gene_scores': factor_scores[top_gene_indices],
                    'program_id': f"{patient_id}_K{k}_F{factor_idx}"
                }
                
                k_programs.append(program)
                all_programs.append(program)
            
            patient_programs[patient_id][k] = k_programs
    
    print(f"✅ Generated {len(all_programs)} total programs from {len(valid_patients)} patients")
    
    # Step 2: Apply robustness criteria
    print("\n🔍 Step 2: Applying robustness criteria...")
    
    # Criterion 1: 70% overlap with different K in same patient
    print("  Applying criterion 1: Cross-K stability within patients...")
    stable_programs = []
    
    for program in all_programs:
        patient_id = program['patient_id']
        k_value = program['k_value']
        genes = program['genes']
        
        # Check overlap with other K values for same patient
        cross_k_overlaps = []
        for other_k in K_RANGE:
            if other_k != k_value and other_k in patient_programs[patient_id]:
                for other_program in patient_programs[patient_id][other_k]:
                    overlap = get_overlap_count(genes, other_program['genes'])
                    cross_k_overlaps.append(overlap)
        
        # Check if any overlap meets the 70% threshold (35/50 genes)
        if any(overlap >= 35 for overlap in cross_k_overlaps):
            stable_programs.append(program)
    
    print(f"    Programs passing criterion 1: {len(stable_programs)}")
    
    # Criterion 2: 20% overlap with programs in other patients
    print("  Applying criterion 2: Cross-patient similarity...")
    cross_patient_programs = []
    
    for program in stable_programs:
        patient_id = program['patient_id']
        genes = program['genes']
        
        # Check overlap with programs from other patients
        cross_patient_overlaps = []
        for other_program in stable_programs:
            if other_program['patient_id'] != patient_id:
                overlap = get_overlap_count(genes, other_program['genes'])
                cross_patient_overlaps.append(overlap)
        
        # Check if any overlap meets the 20% threshold (10/50 genes)
        if any(overlap >= 10 for overlap in cross_patient_overlaps):
            cross_patient_programs.append(program)
    
    print(f"    Programs passing criterion 2: {len(cross_patient_programs)}")
    
    # Criterion 3: <20% overlap with other programs in same patient + selection
    print("  Applying criterion 3: Intra-patient redundancy removal...")
    
    # Rank programs by cross-patient similarity
    program_similarities = []
    for program in cross_patient_programs:
        patient_id = program['patient_id']
        genes = program['genes']
        
        # Calculate average similarity to programs in other patients
        similarities = []
        for other_program in cross_patient_programs:
            if other_program['patient_id'] != patient_id:
                overlap = get_overlap_count(genes, other_program['genes'])
                similarities.append(overlap)
        
        avg_similarity = np.mean(similarities) if similarities else 0
        program_similarities.append((program, avg_similarity))
    
    # Sort by decreasing cross-patient similarity
    program_similarities.sort(key=lambda x: x[1], reverse=True)
    
    # Select programs, avoiding redundancy within patients
    robust_programs = []
    selected_patient_programs = {patient: [] for patient in valid_patients}
    
    for program, similarity in program_similarities:
        patient_id = program['patient_id']
        genes = program['genes']
        
        # Check overlap with already selected programs from same patient
        intra_patient_overlaps = []
        for selected_program in selected_patient_programs[patient_id]:
            overlap = get_overlap_count(genes, selected_program['genes'])
            intra_patient_overlaps.append(overlap)
        
        # Select if no high overlap with selected programs from same patient
        if not any(overlap >= 10 for overlap in intra_patient_overlaps):  # <20% overlap
            robust_programs.append(program)
            selected_patient_programs[patient_id].append(program)
    
    print(f"    Final robust programs: {len(robust_programs)}")
    
    # Step 3: Hierarchical clustering to identify Meta-Programs (MPs)
    print("\n🌳 Step 3: Hierarchical clustering to identify Meta-Programs...")
    
    if len(robust_programs) < 2:
        print("  ⚠️ Too few robust programs for clustering")
        return None
    
    # Create similarity matrix using Jaccard similarity
    n_programs = len(robust_programs)
    similarity_matrix = np.zeros((n_programs, n_programs))
    
    for i in range(n_programs):
        for j in range(n_programs):
            overlap = calculate_gene_overlap(robust_programs[i]['genes'], robust_programs[j]['genes'])
            similarity_matrix[i, j] = overlap
    
    # Convert to distance matrix
    distance_matrix = 1 - similarity_matrix
    
    # Hierarchical clustering
    condensed_distances = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_distances, method='ward')
    
    # Determine optimal number of clusters (Meta-Programs)
    # Try different numbers of clusters and evaluate
    max_clusters = min(MAX_CLUSTER_NUM, len(robust_programs))
    cluster_solutions = {}
    
    for n_clusters in range(2, max_clusters + 1):
        clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
        
        # Calculate within-cluster similarity
        within_cluster_similarities = []
        for cluster_id in range(1, n_clusters + 1):
            cluster_programs = [i for i, c in enumerate(clusters) if c == cluster_id]
            if len(cluster_programs) > 1:
                cluster_similarities = []
                for i, j in combinations(cluster_programs, 2):
                    cluster_similarities.append(similarity_matrix[i, j])
                within_cluster_similarities.extend(cluster_similarities)
        
        avg_within_similarity = np.mean(within_cluster_similarities) if within_cluster_similarities else 0
        cluster_solutions[n_clusters] = {
            'clusters': clusters,
            'within_similarity': avg_within_similarity,
            'n_programs_per_cluster': [np.sum(clusters == i) for i in range(1, n_clusters + 1)]
        }
    
    # Select optimal clustering (highest within-cluster similarity)
    optimal_n_clusters = max(cluster_solutions.keys(), 
                           key=lambda k: cluster_solutions[k]['within_similarity'])
    optimal_clusters = cluster_solutions[optimal_n_clusters]['clusters']
    
    print(f"    Optimal number of Meta-Programs: {optimal_n_clusters}")
    print(f"    Programs per MP: {cluster_solutions[optimal_n_clusters]['n_programs_per_cluster']}")
    
    # Step 4: Generate Meta-Program signatures
    print("\n🧬 Step 4: Generating Meta-Program signatures...")
    
    meta_programs = {}
    
    for mp_id in range(1, optimal_n_clusters + 1):
        mp_program_indices = [i for i, c in enumerate(optimal_clusters) if c == mp_id]
        mp_programs = [robust_programs[i] for i in mp_program_indices]
        
        print(f"  Meta-Program {mp_id}: {len(mp_programs)} constituent programs")
        
        # Collect all genes from constituent programs
        all_mp_genes = []
        for program in mp_programs:
            all_mp_genes.extend(program['genes'])
        
        # Count gene occurrences
        gene_counts = Counter(all_mp_genes)
        
        # Select genes occurring in ≥25% of constituent programs
        min_occurrence = max(1, len(mp_programs) * 0.25)
        candidate_genes = [gene for gene, count in gene_counts.items() if count >= min_occurrence]
        
        # Calculate gene scores (average across programs where gene appears)
        gene_scores = {}
        for gene in candidate_genes:
            scores = []
            for program in mp_programs:
                if gene in program['genes']:
                    gene_idx = program['genes'].index(gene)
                    scores.append(program['gene_scores'][gene_idx])
            gene_scores[gene] = np.mean(scores)
        
        # Sort by score and take top genes
        sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
        
        # Limit to reasonable signature size
        max_signature_size = min(50, len(sorted_genes))
        signature_genes = [gene for gene, score in sorted_genes[:max_signature_size]]
        signature_scores = [score for gene, score in sorted_genes[:max_signature_size]]
        
        meta_programs[mp_id] = {
            'mp_id': mp_id,
            'constituent_programs': mp_programs,
            'signature_genes': signature_genes,
            'signature_scores': signature_scores,
            'n_programs': len(mp_programs),
            'patients_represented': list(set(p['patient_id'] for p in mp_programs))
        }
        
        print(f"    MP{mp_id} signature: {len(signature_genes)} genes, "
              f"{len(meta_programs[mp_id]['patients_represented'])} samples")
    
    # Step 5: Calculate Meta-Program scores for all cells
    print("\n📊 Step 5: Calculating Meta-Program scores for all cells...")
    
    # Calculate MP scores using scanpy
    for mp_id, mp_data in meta_programs.items():
        signature_genes = mp_data['signature_genes']
        available_signature_genes = [g for g in signature_genes if g in malignant_cells.var_names]
        
        if len(available_signature_genes) >= 5:  # Minimum genes for scoring
            sc.tl.score_genes(malignant_cells, available_signature_genes, 
                             score_name=f'MP{mp_id}_score', use_raw=False)
            print(f"    MP{mp_id}: {len(available_signature_genes)} genes used for scoring")
    
    # Step 6: Assign cells to Meta-Programs
    print("\n🎯 Step 6: Assigning cells to Meta-Programs...")
    
    # Get MP score columns
    mp_score_cols = [f'MP{mp_id}_score' for mp_id in meta_programs.keys() 
                     if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
    if len(mp_score_cols) >= 2:
        # Normalize scores by mean subtraction
        for col in mp_score_cols:
            malignant_cells.obs[f'{col}_normalized'] = (malignant_cells.obs[col] - 
                                                       malignant_cells.obs[col].mean())
        
        # Assign cells to MPs
        normalized_cols = [f'{col}_normalized' for col in mp_score_cols]
        mp_assignments = []
        
        for idx in range(malignant_cells.n_obs):
            scores = [malignant_cells.obs.iloc[idx][col] for col in normalized_cols]
            
            if len(scores) >= 2:
                max_score_idx = np.argmax(scores)
                max_score = scores[max_score_idx]
                
                # Get second highest score
                scores_sorted = sorted(scores, reverse=True)
                second_highest = scores_sorted[1] if len(scores_sorted) > 1 else 0
                
                # Assign if highest score > 85% threshold of second highest
                if second_highest == 0 or max_score > 0.85 * abs(second_highest):
                    mp_assignments.append(f'MP{max_score_idx + 1}')
                else:
                    mp_assignments.append('Unresolved')
            else:
                mp_assignments.append('Unresolved')
        
        malignant_cells.obs['MP_assignment'] = mp_assignments
        malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].astype('category')
        
        # Print assignment statistics
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        print("    MP assignment distribution:")
        for mp, count in assignment_counts.items():
            percentage = count / len(malignant_cells) * 100
            print(f"      {mp}: {count} cells ({percentage:.1f}%)")
    
    # Save results
    print("\n💾 Saving robust NMF results...")
    
    robust_nmf_results = {
        'all_programs': all_programs,
        'robust_programs': robust_programs,
        'meta_programs': meta_programs,
        'clustering_info': {
            'optimal_n_clusters': optimal_n_clusters,
            'similarity_matrix': similarity_matrix,
            'linkage_matrix': linkage_matrix,
            'cluster_assignments': optimal_clusters
        },
        'parameters': {
            'K_range': list(K_RANGE),
            'top_genes_per_program': TOP_GENES_PER_PROGRAM,
            'min_cells_per_patient': MIN_CELLS_PER_PATIENT,
            'same_patient_overlap_threshold': SAME_PATIENT_OVERLAP_THRESHOLD,
            'cross_patient_overlap_threshold': CROSS_PATIENT_OVERLAP_THRESHOLD,
            'max_intra_patient_overlap': MAX_INTRA_PATIENT_OVERLAP
        }
    }
    
    with open(os.path.join(save_path, "robust_metabolism_nmf_results.pkl"), 'wb') as f:
        pickle.dump(robust_nmf_results, f)
    
    print(f"✅ Robust NMF analysis completed!")
    print(f"📊 Summary:")
    print(f"   - {len(all_programs)} total programs generated")
    print(f"   - {len(robust_programs)} robust programs identified")  
    print(f"   - {optimal_n_clusters} Meta-Programs discovered")
    print(f"   - {len(valid_patients)} patients analyzed")
    
    return robust_nmf_results

# Manual Meta-Program Merger Function
def merge_meta_programs(malignant_cells, robust_results, 
                       source_mp=7, target_mp=1, 
                       save_path=None, update_scoring=True):
    """
    Manually merge one Meta-Program into another
    
    Parameters:
    -----------
    malignant_cells : AnnData
        Annotated data object with cells
    robust_results : dict
        Results from robust_nmf_analysis
    source_mp : int
        MP to be merged (will be eliminated), default 7
    target_mp : int
        MP to merge into (will be expanded), default 1
    save_path : str
        Path to save updated results
    update_scoring : bool
        Whether to recalculate MP scores after merger
    
    Returns:
    --------
    dict : Updated robust_results
    AnnData : Updated malignant_cells
    """
    
    print(f"🔄 Merging MP{source_mp} into MP{target_mp}...")
    print("="*50)
    
    # Create deep copy to avoid modifying original data
    import copy
    updated_results = copy.deepcopy(robust_results)
    
    # =========================================================================
    # STEP 1: Update Meta-Programs Dictionary
    # =========================================================================
    
    print(f"📊 Step 1: Updating meta-programs dictionary...")
    
    meta_programs = updated_results['meta_programs']
    
    if source_mp not in meta_programs or target_mp not in meta_programs:
        print(f"❌ Error: MP{source_mp} or MP{target_mp} not found in meta_programs!")
        available_mps = list(meta_programs.keys())
        print(f"Available MPs: {available_mps}")
        return None, None
    
    # Get source and target MP data
    source_mp_data = meta_programs[source_mp]
    target_mp_data = meta_programs[target_mp]
    
    print(f"   MP{source_mp}: {len(source_mp_data['signature_genes'])} genes, "
          f"{source_mp_data['n_programs']} programs")
    print(f"   MP{target_mp}: {len(target_mp_data['signature_genes'])} genes, "
          f"{target_mp_data['n_programs']} programs")
    
    # Merge constituent programs
    merged_programs = target_mp_data['constituent_programs'] + source_mp_data['constituent_programs']
    
    # Merge signature genes with proper scoring
    merged_gene_data = merge_signature_genes(
        target_mp_data['signature_genes'], target_mp_data['signature_scores'],
        source_mp_data['signature_genes'], source_mp_data['signature_scores']
    )
    
    # Update target MP with merged data
    meta_programs[target_mp] = {
        'mp_id': target_mp,
        'constituent_programs': merged_programs,
        'signature_genes': merged_gene_data['genes'],
        'signature_scores': merged_gene_data['scores'],
        'n_programs': len(merged_programs),
        'patients_represented': list(set(
            target_mp_data['patients_represented'] + 
            source_mp_data['patients_represented']
        )),
        'merged_from': [source_mp],  # Track merger history
        'merge_timestamp': pd.Timestamp.now().isoformat()
    }
    
    # Remove source MP
    del meta_programs[source_mp]
    
    print(f"   ✅ Merged MP{target_mp}: {len(merged_gene_data['genes'])} genes, "
          f"{len(merged_programs)} programs, "
          f"{len(meta_programs[target_mp]['patients_represented'])} samples")
    
    # =========================================================================
    # STEP 2: Update Clustering Information
    # =========================================================================
    
    print(f"🌳 Step 2: Updating clustering information...")
    
    clustering_info = updated_results['clustering_info']
    
    # Update optimal number of clusters
    old_n_clusters = clustering_info['optimal_n_clusters']
    new_n_clusters = old_n_clusters - 1
    clustering_info['optimal_n_clusters'] = new_n_clusters
    
    # Update cluster assignments
    old_assignments = clustering_info['cluster_assignments'].copy()
    new_assignments = []
    
    for assignment in old_assignments:
        if assignment == source_mp:
            new_assignments.append(target_mp)  # Reassign source to target
        elif assignment > source_mp:
            new_assignments.append(assignment - 1)  # Shift down by 1
        else:
            new_assignments.append(assignment)  # Keep unchanged
    
    clustering_info['cluster_assignments'] = np.array(new_assignments)
    
    print(f"   ✅ Updated clusters: {old_n_clusters} → {new_n_clusters}")
    
    # =========================================================================
    # STEP 3: Renumber Meta-Programs (Optional but Recommended)
    # =========================================================================
    
    print(f"🔢 Step 3: Renumbering Meta-Programs...")
    
    # Create mapping for renumbering
    old_mp_ids = sorted([mp_id for mp_id in meta_programs.keys() if mp_id != source_mp])
    new_mp_mapping = {}
    
    new_id = 1
    for old_id in old_mp_ids:
        new_mp_mapping[old_id] = new_id
        new_id += 1
    
    # Renumber meta_programs dictionary
    new_meta_programs = {}
    for old_id, new_id in new_mp_mapping.items():
        mp_data = meta_programs[old_id].copy()
        mp_data['mp_id'] = new_id
        new_meta_programs[new_id] = mp_data
    
    updated_results['meta_programs'] = new_meta_programs
    
    print(f"   ✅ Renumbered MPs: {dict(new_mp_mapping)}")
    
    # =========================================================================
    # STEP 4: Update Cell Annotations
    # =========================================================================
    
    print(f"🧬 Step 4: Updating cell annotations...")
    
    # Update MP scores in malignant_cells
    if update_scoring:
        print("   Recalculating MP scores...")
        
        # Remove old MP scores
        old_score_cols = [col for col in malignant_cells.obs.columns 
                         if col.startswith('MP') and '_score' in col]
        for col in old_score_cols:
            del malignant_cells.obs[col]
        
        # Calculate new MP scores
        for new_id, mp_data in new_meta_programs.items():
            signature_genes = mp_data['signature_genes']
            available_genes = [g for g in signature_genes if g in malignant_cells.var_names]
            
            if len(available_genes) >= 5:
                sc.tl.score_genes(malignant_cells, available_genes, 
                                 score_name=f'MP{new_id}_score', use_raw=False)
                print(f"     MP{new_id}: {len(available_genes)} genes used for scoring")
    
    # Update MP assignments
    if 'MP_assignment' in malignant_cells.obs.columns:
        print("   Updating MP assignments...")
        
        old_assignments = malignant_cells.obs['MP_assignment'].copy()
        new_assignments = []
        
        for assignment in old_assignments:
            if assignment == f'MP{source_mp}':
                new_assignments.append(f'MP{new_mp_mapping[target_mp]}')
            elif assignment.startswith('MP'):
                try:
                    old_mp_num = int(assignment.replace('MP', ''))
                    if old_mp_num in new_mp_mapping:
                        new_assignments.append(f'MP{new_mp_mapping[old_mp_num]}')
                    else:
                        new_assignments.append('Unresolved')
                except:
                    new_assignments.append(assignment)  # Keep as is if can't parse
            else:
                new_assignments.append(assignment)  # Keep non-MP assignments
        
        malignant_cells.obs['MP_assignment'] = new_assignments
        malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].astype('category')
        
        # Print updated assignment statistics
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        print("   Updated MP assignment distribution:")
        for mp, count in assignment_counts.items():
            percentage = count / len(malignant_cells) * 100
            print(f"     {mp}: {count} cells ({percentage:.1f}%)")
    
    # =========================================================================
    # STEP 5: Update Normalized Scores (if they exist)
    # =========================================================================
    
    print(f"📏 Step 5: Updating normalized scores...")
    
    # Remove old normalized scores
    normalized_cols = [col for col in malignant_cells.obs.columns 
                      if '_score_normalized' in col]
    for col in normalized_cols:
        del malignant_cells.obs[col]
    
    # Recalculate normalized scores
    new_score_cols = [f'MP{mp_id}_score' for mp_id in new_meta_programs.keys() 
                     if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
    if len(new_score_cols) >= 2:
        for col in new_score_cols:
            malignant_cells.obs[f'{col}_normalized'] = (
                malignant_cells.obs[col] - malignant_cells.obs[col].mean()
            )
        print(f"   ✅ Recalculated {len(new_score_cols)} normalized scores")
    
    # =========================================================================
    # STEP 6: Add Merger Metadata
    # =========================================================================
    
    print(f"📝 Step 6: Adding merger metadata...")
    
    # Add merger information to results
    if 'merger_history' not in updated_results:
        updated_results['merger_history'] = []
    
    merger_record = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'source_mp': source_mp,
        'target_mp': target_mp,
        'final_target_mp': new_mp_mapping[target_mp],
        'old_n_clusters': old_n_clusters,
        'new_n_clusters': new_n_clusters,
        'mp_mapping': new_mp_mapping,
        'merged_gene_count': len(merged_gene_data['genes']),
        'merged_program_count': len(merged_programs)
    }
    
    updated_results['merger_history'].append(merger_record)
    
    print(f"   ✅ Merger recorded in metadata")
    
    # =========================================================================
    # STEP 7: Save Updated Results
    # =========================================================================
    
    if save_path:
        print(f"💾 Step 7: Saving updated results...")
        
        # Save updated robust results
        with open(os.path.join(save_path, "robust_nmf_results_merged.pkl"), 'wb') as f:
            pickle.dump(updated_results, f)
        
        # Save updated malignant cells
        malignant_cells.write(os.path.join(save_path, "malignant_cells_merged_mps.h5ad"))
        
        # Save merger summary
        merger_summary = pd.DataFrame([merger_record])
        merger_summary.to_csv(os.path.join(save_path, "mp_merger_summary.csv"), index=False)
        
        print(f"   ✅ Results saved to {save_path}")
    
    # =========================================================================
    # STEP 8: Print Summary
    # =========================================================================
    
    print(f"\n🎉 Meta-Program Merger Completed!")
    print("="*50)
    print(f"✅ Merged MP{source_mp} into MP{target_mp}")
    print(f"📊 Meta-Programs: {old_n_clusters} → {new_n_clusters}")
    print(f"🧬 Final MP{new_mp_mapping[target_mp]} signature: {len(merged_gene_data['genes'])} genes")
    print(f"🔗 Merged programs: {len(merged_programs)} total")
    print(f"📍 Samples represented: {len(new_meta_programs[new_mp_mapping[target_mp]]['patients_represented'])}")
    
    return updated_results, malignant_cells

def merge_signature_genes(target_genes, target_scores, source_genes, source_scores, 
                         max_genes=50, score_weight=0.7):
    """
    Merge signature genes from two Meta-Programs with proper scoring
    
    Parameters:
    -----------
    target_genes, source_genes : list
        Gene lists from target and source MPs
    target_scores, source_scores : list  
        Corresponding gene scores
    max_genes : int
        Maximum number of genes in merged signature
    score_weight : float
        Weight for averaging scores (0.5 = equal weight)
    """
    
    # Create gene score dictionaries
    target_gene_scores = dict(zip(target_genes, target_scores))
    source_gene_scores = dict(zip(source_genes, source_scores))
    
    # Merge genes with weighted averaging for overlapping genes
    merged_gene_scores = {}
    
    # Add target genes
    for gene, score in target_gene_scores.items():
        merged_gene_scores[gene] = score
    
    # Add source genes (average scores for overlapping genes)
    for gene, score in source_gene_scores.items():
        if gene in merged_gene_scores:
            # Weighted average for overlapping genes
            existing_score = merged_gene_scores[gene]
            merged_gene_scores[gene] = (score_weight * existing_score + 
                                      (1 - score_weight) * score)
        else:
            # Add new gene from source
            merged_gene_scores[gene] = score
    
    # Sort by score and select top genes
    sorted_genes = sorted(merged_gene_scores.items(), key=lambda x: x[1], reverse=True)
    
    final_genes = [gene for gene, score in sorted_genes[:max_genes]]
    final_scores = [score for gene, score in sorted_genes[:max_genes]]
    
    return {
        'genes': final_genes,
        'scores': final_scores,
        'n_target_genes': len(target_genes),
        'n_source_genes': len(source_genes),
        'n_overlapping': len(set(target_genes) & set(source_genes)),
        'n_final_genes': len(final_genes)
    }

def validate_merger(updated_results, malignant_cells, save_path=None):
    """
    Validate the merger results and create summary statistics
    """
    
    print("🔍 Validating merger results...")
    
    validation_stats = {}
    
    # Check meta-programs
    meta_programs = updated_results['meta_programs']
    validation_stats['n_meta_programs'] = len(meta_programs)
    validation_stats['mp_ids'] = list(meta_programs.keys())
    
    # Check gene signature sizes
    signature_sizes = {mp_id: len(mp_data['signature_genes']) 
                      for mp_id, mp_data in meta_programs.items()}
    validation_stats['signature_sizes'] = signature_sizes
    
    # Check cell assignments
    if 'MP_assignment' in malignant_cells.obs.columns:
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        validation_stats['assignment_counts'] = assignment_counts.to_dict()
        
        # Check for unresolved cells
        unresolved_pct = (assignment_counts.get('Unresolved', 0) / 
                         len(malignant_cells) * 100)
        validation_stats['unresolved_percentage'] = unresolved_pct
    
    # Check MP scores
    mp_score_cols = [col for col in malignant_cells.obs.columns 
                    if col.startswith('MP') and '_score' in col and 'normalized' not in col]
    validation_stats['available_mp_scores'] = mp_score_cols
    
    # Print validation summary
    print("📊 Validation Summary:")
    print(f"   Meta-Programs: {validation_stats['n_meta_programs']}")
    print(f"   MP IDs: {validation_stats['mp_ids']}")
    print(f"   Signature sizes: {signature_sizes}")
    print(f"   Available MP scores: {len(mp_score_cols)}")
    
    if 'assignment_counts' in validation_stats:
        print(f"   Cell assignments:")
        for mp, count in validation_stats['assignment_counts'].items():
            pct = count / len(malignant_cells) * 100
            print(f"     {mp}: {count} ({pct:.1f}%)")
        print(f"   Unresolved: {unresolved_pct:.1f}%")
    
    if save_path:
        validation_df = pd.DataFrame([validation_stats])
        validation_df.to_csv(os.path.join(save_path, "merger_validation.csv"), index=False)
    
    return validation_stats

def split_meta_program(malignant_cells, robust_results, 
                      mp_to_split, n_splits=2, 
                      split_method='expression_clustering',
                      save_path=None, update_scoring=True,
                      min_cells_per_split=50, min_genes_per_split=10):
    """
    Split a Meta-Program into multiple distinct programs
    
    Parameters:
    -----------
    malignant_cells : AnnData
        Annotated data object with cells
    robust_results : dict
        Results from robust_nmf_analysis
    mp_to_split : int
        MP to be split
    n_splits : int
        Number of programs to split into
    split_method : str
        Method for splitting ('expression_clustering', 'program_clustering', 'gene_clustering')
    save_path : str
        Path to save updated results
    update_scoring : bool
        Whether to recalculate MP scores after split
    min_cells_per_split : int
        Minimum cells required per split
    min_genes_per_split : int
        Minimum genes required per split
    
    Returns:
    --------
    dict : Updated robust_results
    AnnData : Updated malignant_cells
    """
    
    print(f"✂️ Splitting MP{mp_to_split} into {n_splits} programs...")
    print("="*60)
    
    # Create deep copy to avoid modifying original data
    import copy
    updated_results = copy.deepcopy(robust_results)
    
    # =========================================================================
    # STEP 1: Validate Input and Extract MP Data
    # =========================================================================
    
    print(f"📊 Step 1: Validating input and extracting MP data...")
    
    meta_programs = updated_results['meta_programs']
    
    if mp_to_split not in meta_programs:
        print(f"❌ Error: MP{mp_to_split} not found in meta_programs!")
        available_mps = list(meta_programs.keys())
        print(f"Available MPs: {available_mps}")
        return None, None
    
    mp_data = meta_programs[mp_to_split]
    
    print(f"   MP{mp_to_split}: {len(mp_data['signature_genes'])} genes, "
          f"{mp_data['n_programs']} programs, "
          f"{len(mp_data['patients_represented'])} samples")
    
    # Check if MP has enough data to split
    if mp_data['n_programs'] < n_splits:
        print(f"❌ Error: MP{mp_to_split} has only {mp_data['n_programs']} programs, "
              f"cannot split into {n_splits}")
        return None, None
    
    # Get cells assigned to this MP
    mp_cells_mask = malignant_cells.obs['MP_assignment'] == f'MP{mp_to_split}'
    n_mp_cells = mp_cells_mask.sum()
    
    print(f"   Cells assigned to MP{mp_to_split}: {n_mp_cells}")
    
    if n_mp_cells < min_cells_per_split * n_splits:
        print(f"❌ Error: Not enough cells for splitting. Need {min_cells_per_split * n_splits}, "
              f"have {n_mp_cells}")
        return None, None
    
    # =========================================================================
    # STEP 2: Perform Splitting Based on Selected Method
    # =========================================================================
    
    print(f"✂️ Step 2: Performing split using '{split_method}' method...")
    
    if split_method == 'program_clustering':
        split_results = split_by_program_clustering(
            mp_data, n_splits, min_genes_per_split
        )
    elif split_method == 'gene_clustering':
        split_results = split_by_gene_clustering(
            malignant_cells, mp_data, mp_cells_mask, n_splits, 
            min_genes_per_split, save_path
        )
    else:
        print(f"❌ Error: Unknown split method '{split_method}'")
        return None, None
    
    if split_results is None:
        print(f"❌ Splitting failed with method '{split_method}'")
        return None, None
    
    # =========================================================================
    # STEP 3: Create New Meta-Programs from Split Results
    # =========================================================================
    
    print(f"🧬 Step 3: Creating new Meta-Programs from split results...")
    
    # Remove original MP
    del meta_programs[mp_to_split]
    
    # Find the next available MP ID
    existing_mp_ids = list(meta_programs.keys())
    next_mp_id = max(existing_mp_ids) + 1 if existing_mp_ids else 1
    
    new_mp_data = {}
    new_mp_mapping = {}  # Track old MP -> new MPs
    
    for i, split_data in enumerate(split_results):
        new_mp_id = next_mp_id + i
        
        # Create new MP data structure
        new_mp_data[new_mp_id] = {
            'mp_id': new_mp_id,
            'constituent_programs': split_data['programs'],
            'signature_genes': split_data['signature_genes'],
            'signature_scores': split_data['signature_scores'],
            'n_programs': len(split_data['programs']),
            'patients_represented': split_data['patients_represented'],
            'split_from': mp_to_split,
            'split_method': split_method,
            'split_index': i,
            'split_timestamp': pd.Timestamp.now().isoformat()
        }
        
        print(f"   MP{new_mp_id}: {len(split_data['signature_genes'])} genes, "
              f"{len(split_data['programs'])} programs, "
              f"{len(split_data['patients_represented'])} samples")
    
    # Add new MPs to meta_programs
    meta_programs.update(new_mp_data)
    new_mp_mapping[mp_to_split] = list(new_mp_data.keys())
    
    # =========================================================================
    # STEP 4: Update Clustering Information
    # =========================================================================
    
    print(f"🌳 Step 4: Updating clustering information...")
    
    clustering_info = updated_results['clustering_info']
    
    # Update number of clusters
    old_n_clusters = clustering_info['optimal_n_clusters']
    new_n_clusters = old_n_clusters + n_splits - 1
    clustering_info['optimal_n_clusters'] = new_n_clusters
    
    # Update cluster assignments (this is complex for splits)
    old_assignments = clustering_info['cluster_assignments'].copy()
    
    # For simplicity, we'll mark split assignments as requiring recalculation
    clustering_info['cluster_assignments_valid'] = False
    clustering_info['requires_reclustering'] = True
    
    print(f"   ✅ Updated clusters: {old_n_clusters} → {new_n_clusters}")
    print(f"   ⚠️ Cluster assignments marked for recalculation")
    
    # =========================================================================
    # STEP 5: Update Cell Annotations
    # =========================================================================
    
    print(f"🧬 Step 5: Updating cell annotations...")
    
    # Update MP scores
    if update_scoring:
        print("   Recalculating MP scores...")
        
        # Remove old MP scores for split MP
        old_score_cols = [col for col in malignant_cells.obs.columns 
                         if col.startswith(f'MP{mp_to_split}_score')]
        for col in old_score_cols:
            del malignant_cells.obs[col]
        
        # Calculate new MP scores for split MPs
        for new_mp_id, mp_data in new_mp_data.items():
            signature_genes = mp_data['signature_genes']
            available_genes = [g for g in signature_genes if g in malignant_cells.var_names]
            
            if len(available_genes) >= min_genes_per_split:
                sc.tl.score_genes(malignant_cells, available_genes, 
                                 score_name=f'MP{new_mp_id}_score', use_raw=False)
                print(f"     MP{new_mp_id}: {len(available_genes)} genes used for scoring")
    
    # Update MP assignments using new scores
    if 'cell_assignments' in split_results[0]:
        print("   Updating MP assignments based on split results...")
        
        # Update assignments for cells that were in the split MP
        mp_cells_indices = np.where(mp_cells_mask)[0]
        
        for i, split_data in enumerate(split_results):
            new_mp_id = list(new_mp_data.keys())[i]
            
            # Get cell indices for this split
            if 'cell_assignments' in split_data:
                split_cell_indices = split_data['cell_assignments']
                
                # Update assignments
                for cell_idx in split_cell_indices:
                    if cell_idx < len(malignant_cells):
                        malignant_cells.obs.iloc[cell_idx, 
                            malignant_cells.obs.columns.get_loc('MP_assignment')] = f'MP{new_mp_id}'
    
    else:
        # Reassign cells based on highest MP scores
        print("   Reassigning cells based on highest scores...")
        reassign_cells_after_split(malignant_cells, mp_to_split, list(new_mp_data.keys()))
    
    # Update assignment category
    malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].astype('category')
    
    # Print updated assignment statistics
    assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
    print("   Updated MP assignment distribution:")
    for mp, count in assignment_counts.items():
        if mp.startswith('MP') and any(str(mp_id) in mp for mp_id in new_mp_data.keys()):
            percentage = count / len(malignant_cells) * 100
            print(f"     {mp}: {count} cells ({percentage:.1f}%)")
    
    # =========================================================================
    # STEP 6: Add Split Metadata
    # =========================================================================
    
    print(f"📝 Step 6: Adding split metadata...")
    
    # Add split information to results
    if 'split_history' not in updated_results:
        updated_results['split_history'] = []
    
    split_record = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'original_mp': mp_to_split,
        'split_method': split_method,
        'n_splits': n_splits,
        'new_mp_ids': list(new_mp_data.keys()),
        'old_n_clusters': old_n_clusters,
        'new_n_clusters': new_n_clusters,
        'split_gene_counts': [len(mp['signature_genes']) for mp in new_mp_data.values()],
        'split_program_counts': [mp['n_programs'] for mp in new_mp_data.values()]
    }
    
    updated_results['split_history'].append(split_record)
    
    print(f"   ✅ Split recorded in metadata")
    
    # =========================================================================
    # STEP 7: Save Updated Results
    # =========================================================================
    
    if save_path:
        print(f"💾 Step 7: Saving updated results...")
        
        # Save updated robust results
        with open(os.path.join(save_path, "robust_nmf_results_split.pkl"), 'wb') as f:
            pickle.dump(updated_results, f)
        
        # Save updated malignant cells
        malignant_cells.write(os.path.join(save_path, "malignant_cells_split_mps.h5ad"))
        
        # Save split summary
        split_summary = pd.DataFrame([split_record])
        split_summary.to_csv(os.path.join(save_path, "mp_split_summary.csv"), index=False)
        
        print(f"   ✅ Results saved to {save_path}")
    
    # =========================================================================
    # STEP 8: Print Summary
    # =========================================================================
    
    print(f"\n🎉 Meta-Program Split Completed!")
    print("="*60)
    print(f"✂️ Split MP{mp_to_split} into {n_splits} new MPs")
    print(f"📊 Meta-Programs: {old_n_clusters} → {new_n_clusters}")
    print(f"🆕 New MP IDs: {list(new_mp_data.keys())}")
    
    for new_mp_id, mp_data in new_mp_data.items():
        print(f"   MP{new_mp_id}: {len(mp_data['signature_genes'])} genes, "
              f"{mp_data['n_programs']} programs")
    
    return updated_results, malignant_cells

def split_by_program_clustering(mp_data, n_splits, min_genes_per_split):
    """
    Split MP based on clustering of constituent programs
    """
    
    print(f"   🔧 Using program-based clustering...")
    
    constituent_programs = mp_data['constituent_programs']
    
    if len(constituent_programs) < n_splits:
        print(f"   ❌ Not enough constituent programs for clustering")
        return None
    
    # Create similarity matrix between programs
    n_programs = len(constituent_programs)
    similarity_matrix = np.zeros((n_programs, n_programs))
    
    for i in range(n_programs):
        for j in range(n_programs):
            genes_i = set(constituent_programs[i]['genes'])
            genes_j = set(constituent_programs[j]['genes'])
            jaccard_sim = len(genes_i & genes_j) / len(genes_i | genes_j) if genes_i | genes_j else 0
            similarity_matrix[i, j] = jaccard_sim
    
    # Hierarchical clustering
    distance_matrix = 1 - similarity_matrix
    condensed_distances = pdist(distance_matrix)
    linkage_matrix = linkage(condensed_distances, method='ward')
    cluster_labels = fcluster(linkage_matrix, n_splits, criterion='maxclust')
    
    # Create split results
    split_results = []
    
    for cluster_id in range(1, n_splits + 1):
        cluster_programs = [constituent_programs[i] for i, label in enumerate(cluster_labels) 
                           if label == cluster_id]
        
        if len(cluster_programs) == 0:
            continue
        
        # Merge genes from cluster programs
        all_genes = []
        all_scores = []
        
        for program in cluster_programs:
            all_genes.extend(program['genes'])
            all_scores.extend(program.get('gene_scores', [1.0] * len(program['genes'])))
        
        # Count gene frequencies and calculate average scores
        gene_counts = Counter(all_genes)
        gene_score_sums = {}
        
        for gene, score in zip(all_genes, all_scores):
            if gene not in gene_score_sums:
                gene_score_sums[gene] = []
            gene_score_sums[gene].append(score)
        
        # Select genes appearing in multiple programs
        cluster_genes = []
        cluster_scores = []
        
        for gene, count in gene_counts.items():
            if count >= max(1, len(cluster_programs) * 0.3):  # 30% threshold
                avg_score = np.mean(gene_score_sums[gene])
                cluster_genes.append(gene)
                cluster_scores.append(avg_score)
        
        # Sort by score and take top genes
        if len(cluster_genes) >= min_genes_per_split:
            gene_score_pairs = list(zip(cluster_genes, cluster_scores))
            gene_score_pairs.sort(key=lambda x: x[1], reverse=True)
            
            final_genes = [gene for gene, score in gene_score_pairs[:40]]
            final_scores = [score for gene, score in gene_score_pairs[:40]]
            
            split_results.append({
                'signature_genes': final_genes,
                'signature_scores': final_scores,
                'programs': cluster_programs,
                'patients_represented': list(set(p['patient_id'] for p in cluster_programs)),
                'n_programs': len(cluster_programs)
            })
    
    return split_results

def split_by_gene_clustering(malignant_cells, mp_data, mp_cells_mask, 
                           n_splits, min_genes_per_split, save_path=None):
    """
    Split MP based on clustering of signature genes
    """
    
    print(f"   🧬 Using gene-based clustering...")
    
    signature_genes = mp_data['signature_genes']
    
    if len(signature_genes) < min_genes_per_split * n_splits:
        print(f"   ❌ Not enough genes for clustering")
        return None
    
    # Get expression of signature genes in MP cells
    mp_cells = malignant_cells[mp_cells_mask]
    available_genes = [g for g in signature_genes if g in mp_cells.var_names]
    
    if len(available_genes) < min_genes_per_split * n_splits:
        print(f"   ❌ Not enough available genes for clustering")
        return None
    
    # Calculate gene-gene correlation matrix
    X_genes = mp_cells[:, available_genes].X.T  # Genes x Cells
    if hasattr(X_genes, 'toarray'):
        X_genes = X_genes.toarray()
    
    # Compute pairwise correlations
    gene_corr_matrix = np.corrcoef(X_genes)
    gene_corr_matrix = np.nan_to_num(gene_corr_matrix)
    
    # Convert correlation to distance
    gene_distance_matrix = 1 - np.abs(gene_corr_matrix)
    
    # Hierarchical clustering of genes
    condensed_distances = pdist(gene_distance_matrix)
    linkage_matrix = linkage(condensed_distances, method='ward')
    gene_cluster_labels = fcluster(linkage_matrix, n_splits, criterion='maxclust')
    
    # Create split results based on gene clusters
    split_results = []
    
    for cluster_id in range(1, n_splits + 1):
        cluster_gene_indices = np.where(gene_cluster_labels == cluster_id)[0]
        cluster_genes = [available_genes[i] for i in cluster_gene_indices]
        
        if len(cluster_genes) < min_genes_per_split:
            continue
        
        # Calculate scores for genes in this cluster
        cluster_scores = [mp_data['signature_scores'][signature_genes.index(gene)] 
                         for gene in cluster_genes if gene in signature_genes]
        
        # Pad scores if needed
        while len(cluster_scores) < len(cluster_genes):
            cluster_scores.append(1.0)
        
        # Distribute constituent programs
        n_programs_for_cluster = len(mp_data['constituent_programs']) // n_splits
        start_idx = (cluster_id - 1) * n_programs_for_cluster
        end_idx = start_idx + n_programs_for_cluster if cluster_id < n_splits else len(mp_data['constituent_programs'])
        cluster_programs = mp_data['constituent_programs'][start_idx:end_idx]
        
        split_results.append({
            'signature_genes': cluster_genes,
            'signature_scores': cluster_scores,
            'programs': cluster_programs,
            'patients_represented': list(set(p['patient_id'] for p in cluster_programs)),
            'gene_cluster_id': cluster_id
        })
    
    return split_results

def reassign_cells_after_split(malignant_cells, original_mp, new_mp_ids):
    """
    Reassign cells originally assigned to split MP based on new MP scores
    Fixed version that handles Categorical columns properly
    """
    
    # Get cells that were assigned to original MP
    original_assignment_mask = malignant_cells.obs['MP_assignment'] == f'MP{original_mp}'
    
    if original_assignment_mask.sum() == 0:
        print(f"   ⚠️ No cells found assigned to MP{original_mp}")
        return
    
    print(f"   📊 Found {original_assignment_mask.sum()} cells assigned to MP{original_mp}")
    
    # Get scores for new MPs
    new_score_cols = [f'MP{mp_id}_score' for mp_id in new_mp_ids 
                     if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
    if len(new_score_cols) == 0:
        print("   ⚠️ No new MP scores available for reassignment")
        return
    
    print(f"   🔄 Using scores from: {new_score_cols}")
    
    # SOLUTION 1: Convert to object type temporarily, then back to categorical
    # This is the most robust approach
    
    # Store original categorical info
    original_categories = malignant_cells.obs['MP_assignment'].cat.categories.tolist()
    
    # Add new MP categories if they don't exist
    new_categories = [f'MP{mp_id}' for mp_id in new_mp_ids]
    all_categories = list(set(original_categories + new_categories))
    
    # Update categories first
    malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].cat.add_categories(
        [cat for cat in new_categories if cat not in original_categories]
    )
    
    # Now reassign cells based on highest score
    reassignment_counts = {f'MP{mp_id}': 0 for mp_id in new_mp_ids}
    
    for cell_idx in np.where(original_assignment_mask)[0]:
        # Get scores for this cell
        scores = []
        for col in new_score_cols:
            score = malignant_cells.obs.iloc[cell_idx][col]
            scores.append(score)
        
        # Find MP with highest score
        best_mp_idx = np.argmax(scores)
        best_mp_id = new_mp_ids[best_mp_idx]
        best_score = scores[best_mp_idx]
        
        # Assign cell to best MP
        malignant_cells.obs.iloc[cell_idx, 
            malignant_cells.obs.columns.get_loc('MP_assignment')] = f'MP{best_mp_id}'
        
        reassignment_counts[f'MP{best_mp_id}'] += 1
    
    # Print reassignment summary
    print(f"   ✅ Reassignment completed:")
    for mp_name, count in reassignment_counts.items():
        if count > 0:
            percentage = count / original_assignment_mask.sum() * 100
            print(f"     {mp_name}: {count} cells ({percentage:.1f}%)")

def create_split_visualization(X_data, cluster_labels, split_results, title, save_path):
    """
    Create visualization of the split results
    """
    
    # PCA for visualization
    from sklearn.decomposition import PCA
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_data)
    
    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: PCA with cluster colors
    scatter = axes[0].scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels, 
                             cmap='tab10', alpha=0.7)
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    axes[0].set_title(f'{title} - Cell Clustering')
    plt.colorbar(scatter, ax=axes[0])
    
    # Plot 2: Split statistics
    split_names = [f'Split {i+1}' for i in range(len(split_results))]
    gene_counts = [len(result['signature_genes']) for result in split_results]
    
    bars = axes[1].bar(split_names, gene_counts, alpha=0.7)
    axes[1].set_ylabel('Number of Signature Genes')
    axes[1].set_title('Split Results - Gene Counts')
    
    # Add value labels on bars
    for bar, count in zip(bars, gene_counts):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    str(count), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f'{save_path}/{title}_split_visualization.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{save_path}/{title}_split_visualization.png', dpi=300, bbox_inches='tight')
    plt.close()

def export_results_for_r_analysis(malignant_cells, robust_results, save_path):
    """
    Export all necessary results for R statistical analysis and visualization
    """
 
    # Create R analysis directory
    r_analysis_dir = os.path.join(save_path, "NMF_results_for_R")
    os.makedirs(r_analysis_dir, exist_ok=True)
    
    # 1. Clinical metadata with MP assignments
    print("📊 Exporting clinical metadata...")
    clinical_data = malignant_cells.obs.copy()
    
    # Add MP scores if available
    mp_score_cols = [col for col in clinical_data.columns if 'MP' in col and 'score' in col]
    
    # Ensure all required columns are present
    required_cols = ['Patient_ID', 'Sample_ID', 'Treatment_Strategy', 'Microsatellite_Status', 
                     'Response', 'Treatment_Stage', 'Gender', 'Age', 'study', 'Tissue_Type']
    
    # Add MP assignment columns
    mp_cols = [col for col in clinical_data.columns if col.startswith('MP') or 'Meta' in col]
    
    export_cols = [col for col in required_cols + mp_cols + mp_score_cols if col in clinical_data.columns]
    
    clinical_export = clinical_data[export_cols].copy()
    clinical_export.to_csv(os.path.join(r_analysis_dir, "clinical_metadata.csv"), index=True)
    
    # 2. Patient-level summary statistics
    print("📈 Creating patient-level summary...")
    patient_summary = []
    
    for patient_id in clinical_data['Patient_ID'].unique():
        patient_cells = clinical_data[clinical_data['Patient_ID'] == patient_id]
        
        # Basic info
        patient_info = {
            'Patient_ID': patient_id,
            'n_cells': len(patient_cells),
            'Treatment_Strategy': patient_cells['Treatment_Strategy'].iloc[0] if 'Treatment_Strategy' in patient_cells.columns else 'Unknown',
            'Microsatellite_Status': patient_cells['Microsatellite_Status'].iloc[0] if 'Microsatellite_Status' in patient_cells.columns else 'Unknown',
            'Response': patient_cells['Response'].iloc[0] if 'Response' in patient_cells.columns else 'Unknown',
            'Treatment_Stage': patient_cells['Treatment_Stage'].iloc[0] if 'Treatment_Stage' in patient_cells.columns else 'Unknown',
            'Gender': patient_cells['Gender'].iloc[0] if 'Gender' in patient_cells.columns else 'Unknown',
            'Age': patient_cells['Age'].iloc[0] if 'Age' in patient_cells.columns else np.nan,
            'study': patient_cells['study'].iloc[0] if 'study' in patient_cells.columns else 'Unknown'
        }
        
        # MP proportions if available
        if 'MP_assignment' in patient_cells.columns:
            mp_counts = patient_cells['MP_assignment'].value_counts()
            total_cells = len(patient_cells)
            
            for mp in mp_counts.index:
                patient_info[f'{mp}_count'] = mp_counts[mp]
                patient_info[f'{mp}_proportion'] = mp_counts[mp] / total_cells
        
        # MP scores if available
        for col in mp_score_cols:
            if col in patient_cells.columns:
                patient_info[f'{col}_mean'] = patient_cells[col].mean()
                patient_info[f'{col}_median'] = patient_cells[col].median()
                patient_info[f'{col}_std'] = patient_cells[col].std()
        
        patient_summary.append(patient_info)
    
    patient_df = pd.DataFrame(patient_summary)
    patient_df.to_csv(os.path.join(r_analysis_dir, "patient_level_summary.csv"), index=False)
    
    # 3. MP signature genes if available
    print("🧬 Exporting MP signature genes...")
    if robust_results and 'meta_programs' in robust_results:
        mp_signatures = robust_results['meta_programs']
        
        # Create signature matrix
        all_genes = set()
        for mp_name, genes in mp_signatures.items():
            all_genes.update(genes)
        
        signature_matrix = pd.DataFrame(0, index=list(all_genes), columns=list(mp_signatures.keys()))
        
        for mp_name, genes in mp_signatures.items():
            for gene in genes:
                signature_matrix.loc[gene, mp_name] = 1
        
        signature_matrix.to_csv(os.path.join(r_analysis_dir, "mp_signature_genes.csv"))
        
        # Export individual signatures as separate files
        for mp_name, genes in mp_signatures.items():
            pd.DataFrame({'gene': genes}).to_csv(
                os.path.join(r_analysis_dir, f"{mp_name}_signature_genes.csv"), index=False)
    
    # 4. Save expression matrix in efficient sparse format (like scanpy standard)
    print("🔬 Exporting expression data in sparse matrix format...")
    
    # Get key genes (MP signatures + highly variable genes)
    key_genes = set()
    
    # Add MP signature genes
    if robust_results and 'meta_programs' in robust_results:
        for program_ in robust_results['meta_programs'].keys():
            key_genes.update(robust_results['meta_programs'][program_]['signature_genes'])
    
    # Add highly variable genes if available
    if hasattr(malignant_cells, 'var') and 'highly_variable' in malignant_cells.var.columns:
        hvg_genes = malignant_cells.var[malignant_cells.var['highly_variable']].index.tolist()
        key_genes.update(hvg_genes[:1000])  # Top 1000 HVGs
    
    # Filter genes that exist in the data
    available_genes = [g for g in key_genes if g in malignant_cells.var.index]
    
    print(f"   Exporting {len(available_genes)} genes in sparse format...")
    
    if available_genes:
        # Subset AnnData object to key genes
        adata_subset = malignant_cells[:, available_genes].copy()
        
        # Convert expression matrix to sparse CSR format and save as .mtx
        # Ensure matrix is in sparse format
        if scipy.sparse.issparse(adata_subset.X):
            matrix = scipy.sparse.csr_matrix(adata_subset.X)
        else:
            matrix = scipy.sparse.csr_matrix(np.array(adata_subset.X))
        
        # Save expression matrix in Matrix Market format (.mtx)
        scipy.io.mmwrite(os.path.join(r_analysis_dir, "expression_matrix.mtx"), matrix)
        
        # Save gene names (rows)
        adata_subset.var_names.to_series().to_csv(
            os.path.join(r_analysis_dir, "gene_names.csv"), 
            header=['gene_name'], index=False
        )
        
        # Save cell barcodes (columns) 
        adata_subset.obs_names.to_series().to_csv(
            os.path.join(r_analysis_dir, "cell_barcodes.csv"), 
            header=['barcode'], index=False
        )
        
        # Save additional gene metadata if available
        if adata_subset.var.shape[1] > 0:
            adata_subset.var.to_csv(os.path.join(r_analysis_dir, "gene_metadata.csv"))
        
        print(f"   ✅ Sparse matrix saved: {matrix.shape[0]} cells × {matrix.shape[1]} genes")
        print(f"   💾 Matrix sparsity: {(1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1])):.1%}")
    
    # 5. NMF analysis parameters and results summary
    print("📋 Exporting analysis parameters...")
    if robust_results:
        analysis_summary = {
            'n_cells': malignant_cells.n_obs,
            'n_genes': malignant_cells.n_vars,
            'n_patients': clinical_data['Patient_ID'].nunique(),
            'n_meta_programs': len(robust_results.get('meta_programs', {})),
        }
        
        # Add parameter information
        if 'parameters' in robust_results:
            analysis_summary.update(robust_results['parameters'])
        
        pd.DataFrame([analysis_summary]).to_csv(
            os.path.join(r_analysis_dir, "analysis_summary.csv"), index=False)
    
    # 6. Treatment response subset (MSI treated cells)
    print("🎯 Exporting MSI treated cells subset...")
    if 'Treatment_Strategy' in clinical_data.columns and 'Microsatellite_Status' in clinical_data.columns:
        msi_treated_mask = (clinical_data['Microsatellite_Status'] == 'MSI') & \
                          (clinical_data['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx']))
        
        msi_treated_data = clinical_data[msi_treated_mask].copy()
        msi_treated_data.to_csv(os.path.join(r_analysis_dir, "msi_treated_cells.csv"), index=True)
        
        # MSI treated patient-level summary
        msi_patient_summary = patient_df[patient_df['Microsatellite_Status'] == 'MSI'].copy()
        msi_patient_summary.to_csv(os.path.join(r_analysis_dir, "msi_treated_patients.csv"), index=False)
    
    # 6. Create data dictionary
    print("📚 Creating data dictionary...")
    data_dict = {
        'clinical_metadata.csv': 'Cell-level clinical metadata with MP assignments and scores',
        'patient_level_summary.csv': 'Patient-level aggregated data with MP proportions and statistics',
        'mp_signature_genes.csv': 'Binary matrix of genes in each meta-program signature',
        'expression_matrix.mtx': 'Sparse expression matrix in Matrix Market format (cells × genes)',
        'gene_names.csv': 'Gene names corresponding to expression matrix rows',
        'cell_barcodes.csv': 'Cell barcodes corresponding to expression matrix columns',
        'gene_metadata.csv': 'Additional gene metadata (if available)',
        'malignant_cells_complete.h5ad': 'Complete AnnData object with all data',
        'analysis_summary.csv': 'Summary of analysis parameters and basic statistics',
        'msi_treated_cells.csv': 'Subset of MSI treated cells for response analysis',
        'msi_treated_patients.csv': 'Patient-level data for MSI treated patients only'
    }
    
    pd.DataFrame(list(data_dict.items()), columns=['File', 'Description']).to_csv(
        os.path.join(r_analysis_dir, "data_dictionary.csv"), index=False)
    
    # NEW EXPORT 1: MP fraction data for stacked barplot
    print("📊 Exporting MP fraction data for stacked barplot...")
    if 'MP_assignment' in malignant_cells.obs.columns:
        sample_mp_fractions = calculate_mp_fractions_per_sample(malignant_cells, 'MP_assignment')
        sample_mp_fractions.to_csv(os.path.join(r_analysis_dir, "sample_mp_fractions.csv"), index=False)
        print(f"   ✅ MP fractions: {sample_mp_fractions.shape}")
    
    # NEW EXPORT 2: MP marker genes with rankings for heatmap
    print("🧬 Exporting MP marker genes for heatmap...")
    if robust_results and 'meta_programs' in robust_results:
        mp_marker_data = []
        for mp_name, program_ in robust_results['meta_programs'].items():
            for i, gene in enumerate(program_["signature_genes"][:50]):  # Top 50 genes per MP
                mp_marker_data.append({
                    'gene': gene,
                    'Meta_Program': mp_name,
                    'rank': i + 1,
                    'weight': 1.0 - (i * 0.02)  # Decreasing importance
                })
        
        if mp_marker_data:
            mp_marker_df = pd.DataFrame(mp_marker_data)
            mp_marker_df.to_csv(os.path.join(r_analysis_dir, "mp_marker_genes_ranked.csv"), index=False)
            print(f"   ✅ MP marker genes: {len(mp_marker_data)} gene-MP pairs")
    
    # NEW EXPORT 3: Expression subset for heatmap (top 20 genes per MP)
    print("🔥 Exporting expression subset for marker gene heatmap...")
    if robust_results and 'meta_programs' in robust_results:
        # Get top 20 genes per MP
        top_genes = set()
        for program_ in robust_results['meta_programs'].values():
            top_genes.update(program_["signature_genes"][:20])
        
        top_genes_list = list(top_genes)
        
        # Check which genes are available in expression matrix
        available_genes = []
        gene_indices = []
        
        for i, gene in enumerate(malignant_cells.var_names):
            if gene in top_genes_list:
                available_genes.append(gene)
                gene_indices.append(i)
        
        if available_genes:
            # Extract expression subset
            expr_subset = malignant_cells.X[:, gene_indices]
            
            # Convert to dense if sparse
            if hasattr(expr_subset, 'toarray'):
                expr_subset = expr_subset.toarray()
            
            # Create DataFrame (genes as rows, cells as columns)
            expr_subset_df = pd.DataFrame(
                expr_subset.T,  # Transpose so genes are rows
                index=available_genes,
                columns=malignant_cells.obs_names
            )
            
            expr_subset_df.to_csv(os.path.join(r_analysis_dir, "expression_matrix_subset.csv"))
            print(f"   ✅ Expression subset: {expr_subset_df.shape} (genes × cells)")
    
    # Update data dictionary
    data_dict.update({
        'sample_mp_fractions.csv': 'MP fractions per sample for stacked barplot',
        'mp_marker_genes_ranked.csv': 'Ranked marker genes per MP for heatmap',
        'expression_matrix_subset.csv': 'Expression data for top MP marker genes'
    })
    
    print(f"\n✅ Enhanced data export completed!")
    print(f"📁 Total files exported: {len(data_dict)}")
    
    return r_analysis_dir
    
# UMAP Visualization of MP Distribution
def create_mp_umap_visualization(malignant_cells, available_metabolism_genes, save_path, figsize=(15, 10)):
    """
    Create UMAP visualization showing MP distribution across cells
    
    Parameters:
    -----------
    malignant_cells : AnnData
        Malignant cells with MP assignments and scores
    save_path : str
        Directory to save figures
    figsize : tuple
        Figure size
    """
    
    print("🗺️ Creating UMAP visualization of Meta-Program distribution...")
    
    # Create a copy for UMAP analysis
    adata_umap = malignant_cells.copy()
    adata_umap = adata_umap[:,available_metabolism_genes]
    
    # Downsample if too many cells (for computational efficiency)
    if adata_umap.n_obs > 30000:
        print(f"  📉 Downsampling from {adata_umap.n_obs} to 30,000 cells for UMAP...")
        sc.pp.subsample(adata_umap, n_obs=30000, random_state=42)
    
    # Preprocessing for UMAP
    print("  🔄 Preprocessing for UMAP...")
    
    # Normalize and log-transform if not already done
    if adata_umap.X.max() > 50:  # Check if data looks like raw counts
        sc.pp.normalize_total(adata_umap, target_sum=1e4)
        sc.pp.log1p(adata_umap)
    
    # Scale data
    sc.pp.scale(adata_umap, max_value=10)
    
    # Principal component analysis
    print("  📊 Computing PCA...")
    sc.tl.pca(adata_umap, svd_solver='arpack', n_comps=50)
    
    # Compute neighborhood graph
    print("  🕸️ Computing neighborhood graph...")
    sc.pp.neighbors(adata_umap, n_neighbors=15, n_pcs=40)
    
    # Compute UMAP
    print("  🗺️ Computing UMAP embedding...")
    sc.tl.umap(adata_umap, random_state=42)
    
    # Create comprehensive UMAP plots
    print("  🎨 Creating UMAP visualizations...")
    
    # Set up the plotting parameters
    sc.settings.set_figure_params(dpi=300, facecolor='white', figsize=figsize)
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(20, 6))
    
    # 1. UMAP colored by MP assignment
    if 'MP_assignment' in adata_umap.obs.columns:
        ax1 = plt.subplot(1, 3, 1)
        sc.pl.umap(adata_umap, color='MP_assignment', 
                  palette='tab20', size=30, alpha=0.8,
                  title='Meta-Program Assignment', 
                  frameon=False, ax=ax1, show=False)
    
    # 2. UMAP colored by Patient ID (to show batch effects)
    if 'Patient_ID' in adata_umap.obs.columns:
        ax2 = plt.subplot(1, 3, 2)
        sc.pl.umap(adata_umap, color='Patient_ID', 
                  palette='plasma',size=20, alpha=0.6,
                  title='Patient ID', 
                  frameon=False, ax=ax2, show=False, legend_loc=None)
    
    # 3. UMAP colored by Treatment Strategy
    if 'Treatment_Strategy' in adata_umap.obs.columns:
        ax3 = plt.subplot(1, 3, 3)
        sc.pl.umap(adata_umap, color='Treatment_Strategy', 
                  palette='Set2', size=30, alpha=0.8,
                  title='Treatment Strategy', 
                  frameon=False, ax=ax3, show=False)
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/mp_umap.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4-6. UMAP colored by top MP scores
    fig = plt.figure(figsize=(20, 10))
    mp_score_cols = [col for col in adata_umap.obs.columns if 'MP' in col and 'score' in col and not 'normalized' in col]

    all_scores = [] # Calculate global min/max for consistent color scaling
    for mp_col in mp_score_cols:
        all_scores.extend(adata_umap.obs[mp_col].values)
    vmin = min(all_scores)
    vmax = max(all_scores)

    for i, mp_col in enumerate(mp_score_cols):
        # Calculate subplot position (2 rows, 5 columns)
        ax = plt.subplot(2, 5, i + 1)
        
        sc.pl.umap(adata_umap, color=mp_col, 
                size=5, alpha=0.8, cmap='rainbow',  # Purple to yellow colormap
                title=f'{mp_col.replace("_", " ").title()}', 
                frameon=False, ax=ax, show=False,
                vmin=vmin, vmax=vmax)  # Consistent color scaling
        
    plt.tight_layout()
    plt.savefig(f"{save_path}/mp_scores_umap.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create individual high-quality UMAP for MP assignment
    if 'MP_assignment' in adata_umap.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Custom color palette for MPs
        mp_categories = adata_umap.obs['MP_assignment'].cat.categories if hasattr(adata_umap.obs['MP_assignment'], 'cat') else adata_umap.obs['MP_assignment'].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(mp_categories)))
        mp_colors = dict(zip(mp_categories, colors))
        
        sc.pl.umap(adata_umap, color='MP_assignment', 
                  palette=mp_colors, size=40, alpha=0.8,
                  title='Meta-Program Assignment Distribution', 
                  frameon=False, ax=ax, show=False,
                  legend_loc='right margin')
        
        plt.savefig(f"{save_path}/mp_assignment_umap_high_quality.pdf", dpi=300, bbox_inches='tight')
        plt.show()
    
    # Create density plots for each MP
    if 'MP_assignment' in adata_umap.obs.columns:
        mp_categories = adata_umap.obs['MP_assignment'].unique()
        n_mps = len(mp_categories)
        
        fig, axes = plt.subplots(2, (n_mps + 1) // 2, figsize=(5 * ((n_mps + 1) // 2), 10))
        axes = axes.flatten() if n_mps > 2 else [axes] if n_mps == 1 else axes
        
        for i, mp in enumerate(mp_categories):
            if mp != 'Unresolved' and mp != '':  # Skip unresolved cells
                # Create binary annotation for this MP
                adata_umap.obs[f'is_{mp}'] = (adata_umap.obs['MP_assignment'] == mp).astype(int)
                
                sc.pl.umap(adata_umap, color=f'is_{mp}', 
                          size=30, alpha=0.8, cmap='Reds',
                          title=f'{mp} Distribution', 
                          frameon=False, ax=axes[i], show=False)
        
        plt.tight_layout()
        plt.savefig(f"{save_path}/mp_individual_density_umap.pdf", dpi=300, bbox_inches='tight')
        plt.show()
    
    print("  ✅ UMAP visualizations saved!")
    
    # Save UMAP coordinates for potential use in R
    umap_coords = pd.DataFrame({
        'cell_barcode': adata_umap.obs.index,
        'UMAP_1': adata_umap.obsm['X_umap'][:, 0],
        'UMAP_2': adata_umap.obsm['X_umap'][:, 1],
        'MP_assignment': adata_umap.obs['MP_assignment'] if 'MP_assignment' in adata_umap.obs.columns else 'Unknown',
        'Patient_ID': adata_umap.obs['Patient_ID'] if 'Patient_ID' in adata_umap.obs.columns else 'Unknown'
    })
    
    # Add MP scores to coordinates
    for col in mp_score_cols:
        if col in adata_umap.obs.columns:
            umap_coords[col] = adata_umap.obs[col].values
    
    umap_coords.to_csv(f"{save_path}/umap_coordinates_with_mp_data.csv", index=False)
    print(f"  💾 UMAP coordinates saved to: {save_path}/umap_coordinates_with_mp_data.csv")
    
    return adata_umap

# Update the analyze_cluster_distribution function to work with MPs
def analyze_mp_distribution(adata, groupby_var, title_prefix, figurePath):
    """Analyze Meta-Program distribution across a grouping variable"""
    
    if 'MP_assignment' not in adata.obs.columns:
        print(f"⚠️ MP_assignment not found in data")
        return None
    
    # Create contingency table
    ct = pd.crosstab(adata.obs['MP_assignment'], adata.obs[groupby_var])
    ct_normalized = pd.crosstab(adata.obs['MP_assignment'], adata.obs[groupby_var], 
                                normalize='columns') * 100
    
    # Statistical test (Chi-square)
    from scipy.stats import chi2_contingency
    chi2, p_value, dof, expected = chi2_contingency(ct)
    
    print(f"\n{title_prefix}:")
    print(f"Chi-square test: χ² = {chi2:.3f}, p-value = {p_value:.3e}")
    print("Normalized distribution (%):")
    print(ct_normalized.round(1))
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Stacked bar plot
    ct_normalized.T.plot(kind='bar', stacked=True, ax=ax1, 
                        colormap='Set3', alpha=0.8, width=0.8)
    ax1.set_title(f'{title_prefix}', fontweight='bold', fontsize=14)
    ax1.set_xlabel(groupby_var, fontweight='bold')
    ax1.set_ylabel('Percentage', fontweight='bold')
    ax1.legend(title='Meta-Program', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    
    # Heatmap
    sns.heatmap(ct_normalized, annot=True, fmt='.1f', cmap='Blues', ax=ax2,
                cbar_kws={'shrink': 0.8})
    ax2.set_title(f'{title_prefix} - Heatmap', fontweight='bold', fontsize=14)
    ax2.set_xlabel(groupby_var, fontweight='bold')
    ax2.set_ylabel('Meta-Program', fontweight='bold')
    plt.setp(ax2.get_xticklabels(), rotation=45)
    
    plt.tight_layout()
    filename = f"{title_prefix.lower().replace(' ', '_')}_mp_distribution"
    plt.savefig(f"{figurePath}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    return {'contingency_table': ct, 'normalized_table': ct_normalized, 
            'chi2': chi2, 'p_value': p_value}

# Enhanced patient-level analysis for Meta-Programs
def comprehensive_mp_resistance_analysis(adata, save_path):
    """
    Comprehensive analysis to identify therapy-resistant Meta-Programs
    
    Analysis strategy:
    1. Pre-treatment baseline differences (non-pCR vs pCR) - INTRINSIC resistance
    2. Post-treatment changes in non-pCR patients (Post vs Pre) - ACQUIRED resistance  
    3. Differential treatment response (non-pCR Post-Pre vs pCR Post-Pre) - THERAPY-SPECIFIC resistance
    """
    
    print("\n🔬 Comprehensive Meta-Program Resistance Analysis")
    print("="*70)
    
    # Analysis 1: Intrinsic Resistance (Pre-treatment baseline)
    print("\n🎯 Analysis 1: Intrinsic Resistance Patterns (Pre-treatment)")
    print("-" * 60)
    
    intrinsic_results = analyze_intrinsic_resistance(adata, base_mask, save_path)
    
    # Analysis 2: Acquired Resistance (Post vs Pre in non-pCR)
    print("\n🔥 Analysis 2: Acquired Resistance Patterns (Post vs Pre in non-pCR)")
    print("-" * 60)
    
    acquired_results = analyze_acquired_resistance(adata, base_mask, save_path)
    
    # Analysis 3: Therapy-Specific Resistance (Differential response to treatment)
    print("\n⚡ Analysis 3: Therapy-Specific Resistance (Differential treatment response)")
    print("-" * 60)
    
    therapy_specific_results = analyze_therapy_specific_resistance(adata, base_mask, save_path)
    
    # Comprehensive summary and visualization
    print("\n📋 Creating Comprehensive Resistance Summary")
    print("-" * 60)
    
    create_resistance_summary(intrinsic_results, acquired_results, therapy_specific_results, save_path)
    
    return {
        'intrinsic_resistance': intrinsic_results,
        'acquired_resistance': acquired_results, 
        'therapy_specific_resistance': therapy_specific_results
    }

def analyze_intrinsic_resistance(adata, base_mask, save_path):
    """Analysis 1: Pre-treatment baseline differences (Intrinsic resistance)"""
    
    # Filter for pre-treatment samples
    intrinsic_mask = base_mask & (adata.obs['Treatment_Stage'] == 'Pre')
    intrinsic_data = adata[intrinsic_mask].copy()
    
    print(f"Pre-treatment cohort: {intrinsic_data.n_obs:,} cells")
    print("Response distribution:")
    print(intrinsic_data.obs['Response'].value_counts())
    
    return run_mp_comparison_analysis(
        intrinsic_data, 
        comparison_type='intrinsic',
        group_col='Response',
        group1='non_pCR', 
        group2='pCR',
        save_path=save_path,
        description="Pre-treatment non-pCR vs pCR (Intrinsic Resistance)"
    )

def analyze_acquired_resistance(adata, base_mask, save_path):
    """Analysis 2: Post vs Pre changes in non-pCR patients (Acquired resistance)"""
    
    # Filter for non-pCR patients with both Pre and Post samples
    nonpcr_mask = base_mask & (adata.obs['Response'] == 'non_pCR')
    nonpcr_data = adata[nonpcr_mask].copy()
    
    # Check which patients have both Pre and Post samples
    patient_stage_counts = nonpcr_data.obs.groupby('Patient_ID')['Treatment_Stage'].nunique()
    paired_patients = patient_stage_counts[patient_stage_counts >= 2].index.tolist()
    
    print(f"Non-pCR patients with paired samples: {len(paired_patients)}")
    
    if len(paired_patients) == 0:
        print("⚠️ No paired Pre/Post samples found for non-pCR patients")
        return None
    
    # Filter for paired patients only
    paired_mask = nonpcr_data.obs['Patient_ID'].isin(paired_patients)
    paired_nonpcr_data = nonpcr_data[paired_mask].copy()
    
    print(f"Paired non-pCR analysis cohort: {paired_nonpcr_data.n_obs:,} cells")
    print("Treatment stage distribution:")
    print(paired_nonpcr_data.obs['Treatment_Stage'].value_counts())
    
    return run_mp_comparison_analysis(
        paired_nonpcr_data,
        comparison_type='acquired', 
        group_col='Treatment_Stage',
        group1='Post',
        group2='Pre', 
        save_path=save_path,
        description="Post vs Pre in non-pCR patients (Acquired Resistance)",
        paired_analysis=True,
        patient_col='Patient_ID'
    )

def analyze_therapy_specific_resistance(adata, base_mask, save_path):
    """Analysis 3: Differential treatment response (Therapy-specific resistance)"""
    
    # Calculate treatment-induced changes for each patient
    treatment_changes = calculate_treatment_induced_changes(adata, base_mask)
    
    if treatment_changes is None:
        print("⚠️ Insufficient data for therapy-specific resistance analysis")
        return None
    
    print(f"Patients with treatment change data: {len(treatment_changes)}")
    
    # Compare treatment-induced changes between non-pCR vs pCR
    return compare_treatment_changes(treatment_changes, save_path)

def run_mp_comparison_analysis(adata, comparison_type, group_col, group1, group2, 
                              save_path, description, paired_analysis=False, patient_col=None):
    """Generic function to run MP comparison analysis"""
    
    print(f"\n{description}")
    
    # Calculate MP fractions per patient
    patient_mp_data = []
    
    for patient_id in adata.obs['Patient_ID'].unique():
        patient_data = adata[adata.obs['Patient_ID'] == patient_id]
        
        if paired_analysis:
            # For paired analysis, calculate fractions for each treatment stage
            for stage in patient_data.obs[group_col].unique():
                stage_data = patient_data[patient_data.obs[group_col] == stage]
                
                if len(stage_data) < 5:  # Minimum cells per patient-stage
                    continue
                
                mp_counts = stage_data.obs['MP_assignment'].value_counts()
                total_cells = len(stage_data)
                
                for mp in adata.obs['MP_assignment'].unique():
                    if mp != 'Unresolved':
                        fraction = mp_counts.get(mp, 0) / total_cells
                        
                        patient_mp_data.append({
                            'patient_id': patient_id,
                            'mp': mp,
                            'fraction': fraction,
                            'group': stage,
                            'total_cells': total_cells
                        })
        else:
            # For unpaired analysis
            if len(patient_data) < 10:
                continue
                
            group_value = patient_data.obs[group_col].iloc[0]
            mp_counts = patient_data.obs['MP_assignment'].value_counts()
            total_cells = len(patient_data)
            
            for mp in adata.obs['MP_assignment'].unique():
                if mp != 'Unresolved':
                    fraction = mp_counts.get(mp, 0) / total_cells
                    
                    patient_mp_data.append({
                        'patient_id': patient_id,
                        'mp': mp,
                        'fraction': fraction,
                        'group': group_value,
                        'total_cells': total_cells
                    })
    
    if not patient_mp_data:
        print("No valid patient data for analysis")
        return None
    
    patient_df = pd.DataFrame(patient_mp_data)
    
    # Statistical analysis
    from scipy import stats
    mp_results = []
    
    for mp in patient_df['mp'].unique():
        mp_data = patient_df[patient_df['mp'] == mp]
        
        group1_fractions = mp_data[mp_data['group'] == group1]['fraction']
        group2_fractions = mp_data[mp_data['group'] == group2]['fraction']
        
        if len(group1_fractions) >= 3 and len(group2_fractions) >= 3:
            
            if paired_analysis:
                # Paired t-test for paired samples
                # Match patients between groups
                group1_patients = set(mp_data[mp_data['group'] == group1]['patient_id'])
                group2_patients = set(mp_data[mp_data['group'] == group2]['patient_id'])
                common_patients = group1_patients & group2_patients
                
                if len(common_patients) >= 3:
                    paired_group1 = []
                    paired_group2 = []
                    
                    for patient in common_patients:
                        g1_fraction = mp_data[(mp_data['patient_id'] == patient) & 
                                            (mp_data['group'] == group1)]['fraction'].iloc[0]
                        g2_fraction = mp_data[(mp_data['patient_id'] == patient) & 
                                            (mp_data['group'] == group2)]['fraction'].iloc[0]
                        paired_group1.append(g1_fraction)
                        paired_group2.append(g2_fraction)
                    
                    t_stat, p_value = stats.ttest_rel(paired_group1, paired_group2)
                    
                    # Calculate effect size for paired data
                    differences = np.array(paired_group1) - np.array(paired_group2)
                    cohens_d = np.mean(differences) / np.std(differences, ddof=1) if np.std(differences, ddof=1) > 0 else 0
                    
                    mp_results.append({
                        'mp': mp,
                        'group1_mean': np.mean(paired_group1),
                        'group2_mean': np.mean(paired_group2),
                        'group1_std': np.std(paired_group1),
                        'group2_std': np.std(paired_group2),
                        'n_paired_patients': len(common_patients),
                        't_statistic': t_stat,
                        'p_value': p_value,
                        'cohens_d': cohens_d,
                        'test_type': 'paired_t_test'
                    })
                else:
                    continue
            else:
                # Unpaired t-test
                t_stat, p_value = stats.ttest_ind(group1_fractions, group2_fractions)
                
                # Calculate effect size (Cohen's d)
                pooled_std = np.sqrt(((len(group1_fractions) - 1) * np.var(group1_fractions, ddof=1) + 
                                    (len(group2_fractions) - 1) * np.var(group2_fractions, ddof=1)) / 
                                   (len(group1_fractions) + len(group2_fractions) - 2))
                
                cohens_d = (np.mean(group1_fractions) - np.mean(group2_fractions)) / pooled_std if pooled_std > 0 else 0
                
                mp_results.append({
                    'mp': mp,
                    'group1_mean': np.mean(group1_fractions),
                    'group2_mean': np.mean(group2_fractions),
                    'group1_std': np.std(group1_fractions),
                    'group2_std': np.std(group2_fractions),
                    'n_group1_patients': len(group1_fractions),
                    'n_group2_patients': len(group2_fractions),
                    't_statistic': t_stat,
                    'p_value': p_value,
                    'cohens_d': cohens_d,
                    'test_type': 'unpaired_t_test'
                })
            
            # Print results
            print(f"\n{mp}:")
            print(f"  {group1} mean: {np.mean(group1_fractions):.3f} ± {np.std(group1_fractions):.3f}")
            print(f"  {group2} mean: {np.mean(group2_fractions):.3f} ± {np.std(group2_fractions):.3f}")
            print(f"  t-test: t={t_stat:.3f}, p={p_value:.3e}")
            print(f"  Cohen's d: {cohens_d:.3f}")
            
            if p_value < 0.05:
                direction = "higher" if cohens_d > 0 else "lower"
                print(f"  🌟 SIGNIFICANT: {group1} {direction} than {group2}")
                if abs(cohens_d) > 0.5:
                    print(f"  📈 LARGE EFFECT SIZE")
    
    if mp_results:
        mp_results_df = pd.DataFrame(mp_results)
        
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        _, p_corrected, _, _ = multipletests(mp_results_df['p_value'], method='fdr_bh')
        mp_results_df['p_corrected'] = p_corrected
        
        # Add analysis metadata
        mp_results_df['comparison_type'] = comparison_type
        mp_results_df['group1'] = group1
        mp_results_df['group2'] = group2
        
        # Identify significant MPs
        significant_mps = mp_results_df[
            (mp_results_df['p_corrected'] < 0.05) & 
            (abs(mp_results_df['cohens_d']) > 0.5)
        ]
        
        if len(significant_mps) > 0:
            print(f"\n🎉 SIGNIFICANT META-PROGRAMS ({comparison_type}):")
            for _, row in significant_mps.iterrows():
                direction = "higher" if row['cohens_d'] > 0 else "lower"
                print(f"  {row['mp']}: {direction} in {group1}")
                print(f"    {group1}: {row['group1_mean']:.3f} vs {group2}: {row['group2_mean']:.3f}")
                print(f"    Corrected p-value: {row['p_corrected']:.3e}")
                print(f"    Effect size: {row['cohens_d']:.3f}")
        
        # Save results
        filename = f"mp_{comparison_type}_resistance_analysis.csv"
        mp_results_df.to_csv(f"{save_path}/{filename}", index=False)
        
        return {
            'results_df': mp_results_df,
            'significant_mps': significant_mps,
            'comparison_type': comparison_type,
            'description': description
        }
    
    return None

def calculate_treatment_induced_changes(adata, base_mask):
    """Calculate treatment-induced MP changes for each patient"""
    
    # Filter data
    analysis_data = adata[base_mask].copy()
    
    # Find patients with both Pre and Post samples
    patient_stages = analysis_data.obs.groupby('Patient_ID')['Treatment_Stage'].unique()
    paired_patients = []
    
    for patient_id, stages in patient_stages.items():
        if 'Pre' in stages and 'Post' in stages:
            paired_patients.append(patient_id)
    
    if len(paired_patients) == 0:
        return None
    
    treatment_changes = []
    
    for patient_id in paired_patients:
        patient_data = analysis_data[analysis_data.obs['Patient_ID'] == patient_id]
        response = patient_data.obs['Response'].iloc[0]
        
        pre_data = patient_data[patient_data.obs['Treatment_Stage'] == 'Pre']
        post_data = patient_data[patient_data.obs['Treatment_Stage'] == 'Post']
        
        if len(pre_data) >= 5 and len(post_data) >= 5:
            # Calculate MP fractions for Pre and Post
            for mp in analysis_data.obs['MP_assignment'].unique():
                if mp != 'Unresolved':
                    pre_fraction = (pre_data.obs['MP_assignment'] == mp).mean()
                    post_fraction = (post_data.obs['MP_assignment'] == mp).mean()
                    
                    change = post_fraction - pre_fraction
                    
                    treatment_changes.append({
                        'patient_id': patient_id,
                        'mp': mp,
                        'pre_fraction': pre_fraction,
                        'post_fraction': post_fraction,
                        'change': change,
                        'response': response
                    })
    
    return pd.DataFrame(treatment_changes) if treatment_changes else None

def compare_treatment_changes(treatment_changes_df, save_path):
    """Compare treatment-induced changes between non-pCR vs pCR patients"""
    
    from scipy import stats
    change_results = []
    
    for mp in treatment_changes_df['mp'].unique():
        mp_data = treatment_changes_df[treatment_changes_df['mp'] == mp]
        
        nonpcr_changes = mp_data[mp_data['response'] == 'non_pCR']['change']
        pcr_changes = mp_data[mp_data['response'] == 'pCR']['change']
        
        if len(nonpcr_changes) >= 3 and len(pcr_changes) >= 3:
            t_stat, p_value = stats.ttest_ind(nonpcr_changes, pcr_changes)
            
            # Calculate effect size
            pooled_std = np.sqrt(((len(nonpcr_changes) - 1) * np.var(nonpcr_changes, ddof=1) + 
                                (len(pcr_changes) - 1) * np.var(pcr_changes, ddof=1)) / 
                               (len(nonpcr_changes) + len(pcr_changes) - 2))
            
            cohens_d = (np.mean(nonpcr_changes) - np.mean(pcr_changes)) / pooled_std if pooled_std > 0 else 0
            
            change_results.append({
                'mp': mp,
                'nonpcr_change_mean': np.mean(nonpcr_changes),
                'pcr_change_mean': np.mean(pcr_changes),
                'nonpcr_change_std': np.std(nonpcr_changes),
                'pcr_change_std': np.std(pcr_changes),
                'n_nonpcr_patients': len(nonpcr_changes),
                'n_pcr_patients': len(pcr_changes),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d
            })
            
            print(f"\n{mp} (Treatment-induced change):")
            print(f"  Non-pCR change: {np.mean(nonpcr_changes):.3f} ± {np.std(nonpcr_changes):.3f}")
            print(f"  pCR change: {np.mean(pcr_changes):.3f} ± {np.std(pcr_changes):.3f}")
            print(f"  Difference: {np.mean(nonpcr_changes) - np.mean(pcr_changes):.3f}")
            print(f"  t-test: t={t_stat:.3f}, p={p_value:.3e}")
            print(f"  Cohen's d: {cohens_d:.3f}")
    
    if change_results:
        change_results_df = pd.DataFrame(change_results)
        
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        _, p_corrected, _, _ = multipletests(change_results_df['p_value'], method='fdr_bh')
        change_results_df['p_corrected'] = p_corrected
        
        # Save results
        change_results_df.to_csv(f"{save_path}/mp_therapy_specific_resistance.csv", index=False)
        
        # Identify significant therapy-resistant MPs
        therapy_resistant_mps = change_results_df[
            (change_results_df['p_corrected'] < 0.05) & 
            (change_results_df['cohens_d'] > 0.5) &  # More upregulated in non-pCR
            (change_results_df['nonpcr_change_mean'] > 0)  # Actually increased post-treatment
        ]
        
        if len(therapy_resistant_mps) > 0:
            print(f"\n🔥 THERAPY-RESISTANT META-PROGRAMS:")
            for _, row in therapy_resistant_mps.iterrows():
                print(f"  {row['mp']}: Upregulated in non-pCR after treatment")
                print(f"    Non-pCR change: +{row['nonpcr_change_mean']:.3f}")
                print(f"    pCR change: {row['pcr_change_mean']:.3f}")
                print(f"    Corrected p-value: {row['p_corrected']:.3e}")
        
        return {
            'results_df': change_results_df,
            'therapy_resistant_mps': therapy_resistant_mps,
            'comparison_type': 'therapy_specific'
        }
    
    return None

def create_resistance_summary(intrinsic_results, acquired_results, therapy_specific_results, save_path):
    """Create comprehensive resistance analysis summary"""
    
    print("\n📊 Creating Resistance Analysis Summary")
    
    # Combine all significant results
    all_significant_mps = []
    
    if intrinsic_results and 'significant_mps' in intrinsic_results:
        for _, row in intrinsic_results['significant_mps'].iterrows():
            all_significant_mps.append({
                'mp': row['mp'],
                'resistance_type': 'Intrinsic',
                'effect_size': row['cohens_d'],
                'p_corrected': row['p_corrected'],
                'description': f"Enriched in pre-treatment non-pCR"
            })
    
    if acquired_results and 'significant_mps' in acquired_results:
        for _, row in acquired_results['significant_mps'].iterrows():
            all_significant_mps.append({
                'mp': row['mp'],
                'resistance_type': 'Acquired',
                'effect_size': row['cohens_d'],
                'p_corrected': row['p_corrected'],
                'description': f"Upregulated post-treatment in non-pCR"
            })
    
    if therapy_specific_results and 'therapy_resistant_mps' in therapy_specific_results:
        for _, row in therapy_specific_results['therapy_resistant_mps'].iterrows():
            all_significant_mps.append({
                'mp': row['mp'],
                'resistance_type': 'Therapy-Specific',
                'effect_size': row['cohens_d'],
                'p_corrected': row['p_corrected'],
                'description': f"Differentially upregulated in non-pCR vs pCR"
            })
    
    if all_significant_mps:
        summary_df = pd.DataFrame(all_significant_mps)
        summary_df.to_csv(f"{save_path}/comprehensive_resistance_summary.csv", index=False)
        
        print(f"💡 Summary of Resistant Meta-Programs:")
        for resistance_type in summary_df['resistance_type'].unique():
            type_mps = summary_df[summary_df['resistance_type'] == resistance_type]
            print(f"  {resistance_type}: {len(type_mps)} MP(s)")
            for _, row in type_mps.iterrows():
                print(f"    - {row['mp']} (d={row['effect_size']:.3f}, p={row['p_corrected']:.3e})")
    
    print(f"\n✅ Comprehensive resistance analysis completed!")

def comprehensive_msi_vs_mss_resistance_analysis(adata, save_path):
    """
    Comprehensive resistance analysis comparing MSI vs MSS patterns
    
    Analysis strategy:
    1. MSI resistance patterns (previous analysis)
    2. MSS resistance patterns (new analysis)
    3. MSI vs MSS direct comparison
    4. Universal vs specific resistance mechanisms
    """
    
    print("\n🔬 Comprehensive MSI vs MSS Resistance Analysis")
    print("="*70)
    
    # Data overview across microsatellite status
    print("📊 Complete Data Overview:")
    total_treated_mask = adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx'])
    total_treated_data = adata[total_treated_mask].copy()
    
    print(f"Total treated cells: {total_treated_data.n_obs:,}")
    
    # Cross-tabulation of key variables
    print("\nMicrosatellite Status × Treatment Stage × Response:")
    detailed_crosstab = pd.crosstab(
        [total_treated_data.obs['Microsatellite_Status'], total_treated_data.obs['Treatment_Stage']], 
        total_treated_data.obs['Response'], 
        margins=True
    )
    print(detailed_crosstab)
    
    # Check MSS response rates
    mss_data = total_treated_data[total_treated_data.obs['Microsatellite_Status'] == 'MSS']
    if len(mss_data) > 0:
        print(f"\nMSS Response Analysis:")
        mss_response_counts = mss_data.obs['Response'].value_counts()
        print(mss_response_counts)
        
        if 'pCR' in mss_response_counts and 'non_pCR' in mss_response_counts:
            pCR_rate_mss = mss_response_counts['pCR'] / mss_response_counts.sum() * 100
            print(f"MSS pCR rate: {pCR_rate_mss:.1f}%")
        else:
            print("⚠️ MSS response data incomplete")
    
    # Analysis 1: MSI Resistance (Reference - already done)
    print("\n🎯 Analysis 1: MSI Resistance Patterns (Reference)")
    print("-" * 60)
    msi_results = analyze_microsatellite_specific_resistance(
        adata, 'MSI', save_path, "MSI Tumors"
    )
    
    # Analysis 2: MSS Resistance  
    print("\n🔍 Analysis 2: MSS Resistance Patterns")
    print("-" * 60)
    mss_results = analyze_microsatellite_specific_resistance(
        adata, 'MSS', save_path, "MSS Tumors"
    )
    
    # Analysis 3: Direct MSI vs MSS Comparison
    print("\n⚡ Analysis 3: Direct MSI vs MSS Comparison")
    print("-" * 60)
    direct_comparison_results = analyze_msi_vs_mss_programs(adata, save_path)
    
    # Analysis 4: Universal vs Specific Mechanisms
    print("\n🌐 Analysis 4: Universal vs Microsatellite-Specific Mechanisms")
    print("-" * 60)
    universal_vs_specific = identify_universal_vs_specific_resistance(
        msi_results, mss_results, save_path
    )
    
    # Analysis 5: Enhanced Clinical Insights
    print("\n🏥 Analysis 5: Enhanced Clinical Insights")
    print("-" * 60)
    clinical_insights = generate_clinical_insights(
        msi_results, mss_results, direct_comparison_results, save_path
    )
    
    return {
        'msi_resistance': msi_results,
        'mss_resistance': mss_results, 
        'direct_comparison': direct_comparison_results,
        'universal_vs_specific': universal_vs_specific,
        'clinical_insights': clinical_insights
    }

def analyze_microsatellite_specific_resistance(adata, msi_status, save_path, description):
    """Analyze resistance patterns within specific microsatellite status"""
    
    print(f"\n{description} Analysis")
    
    # Filter for specific microsatellite status
    ms_mask = (adata.obs['Microsatellite_Status'] == msi_status) & \
              (adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx']))
    
    ms_data = adata[ms_mask].copy()
    
    print(f"{msi_status} treated cells: {ms_data.n_obs:,}")
    
    if len(ms_data) == 0:
        print(f"⚠️ No {msi_status} treated cells found")
        return None
    
    # Check response distribution
    response_dist = ms_data.obs['Response'].value_counts()
    print(f"{msi_status} Response distribution:")
    print(response_dist)
    
    # Check if we have sufficient data for analysis
    if 'pCR' not in response_dist or 'non_pCR' not in response_dist:
        print(f"⚠️ Insufficient response data for {msi_status} analysis")
        return None
    
    if response_dist['pCR'] < 100 or response_dist['non_pCR'] < 100:
        print(f"⚠️ Low sample sizes for {msi_status}: pCR={response_dist.get('pCR', 0)}, non_pCR={response_dist.get('non_pCR', 0)}")
    
    # Run the three types of resistance analysis
    results = {}
    
    # 1. Intrinsic resistance (Pre-treatment)
    intrinsic_results = analyze_intrinsic_resistance_specific(ms_data, msi_status, save_path)
    results['intrinsic'] = intrinsic_results
    
    # 2. Acquired resistance (Post vs Pre in non-pCR) 
    acquired_results = analyze_acquired_resistance_specific(ms_data, msi_status, save_path)
    results['acquired'] = acquired_results
    
    # 3. Therapy-specific resistance
    therapy_results = analyze_therapy_specific_resistance_specific(ms_data, msi_status, save_path)
    results['therapy_specific'] = therapy_results
    
    return results

def analyze_intrinsic_resistance_specific(ms_data, msi_status, save_path):
    """Intrinsic resistance analysis for specific microsatellite status"""
    
    # Filter for pre-treatment
    pre_mask = ms_data.obs['Treatment_Stage'] == 'Pre'
    pre_data = ms_data[pre_mask].copy()
    
    print(f"\n{msi_status} Pre-treatment analysis:")
    print(f"  Cells: {pre_data.n_obs:,}")
    print(f"  Response distribution: {pre_data.obs['Response'].value_counts().to_dict()}")
    
    if len(pre_data) < 100:
        print(f"  ⚠️ Insufficient pre-treatment {msi_status} data")
        return None
    
    return run_mp_comparison_analysis(
        pre_data,
        comparison_type=f'{msi_status}_intrinsic',
        group_col='Response',
        group1='non_pCR',
        group2='pCR', 
        save_path=save_path,
        description=f"{msi_status} Pre-treatment non-pCR vs pCR",
        paired_analysis=False
    )

def analyze_acquired_resistance_specific(ms_data, msi_status, save_path):
    """Acquired resistance analysis for specific microsatellite status"""
    
    # Filter for non-pCR with paired samples
    nonpcr_data = ms_data[ms_data.obs['Response'] == 'non_pCR'].copy()
    
    # Find patients with both Pre and Post
    patient_stages = nonpcr_data.obs.groupby('Patient_ID')['Treatment_Stage'].unique()
    paired_patients = [pid for pid, stages in patient_stages.items() 
                      if 'Pre' in stages and 'Post' in stages]
    
    print(f"\n{msi_status} Acquired resistance analysis:")
    print(f"  Non-pCR patients with paired samples: {len(paired_patients)}")
    
    if len(paired_patients) < 3:
        print(f"  ⚠️ Insufficient paired {msi_status} data")
        return None
    
    paired_data = nonpcr_data[nonpcr_data.obs['Patient_ID'].isin(paired_patients)].copy()
    
    return run_mp_comparison_analysis(
        paired_data,
        comparison_type=f'{msi_status}_acquired',
        group_col='Treatment_Stage',
        group1='Post',
        group2='Pre',
        save_path=save_path,
        description=f"{msi_status} Post vs Pre in non-pCR",
        paired_analysis=True,
        patient_col='Patient_ID'
    )

def analyze_therapy_specific_resistance_specific(ms_data, msi_status, save_path):
    """Therapy-specific resistance for specific microsatellite status"""
    
    # Calculate treatment-induced changes
    ms_changes = calculate_treatment_induced_changes_ms(ms_data, msi_status)
    
    if ms_changes is None:
        print(f"  ⚠️ No treatment change data for {msi_status}")
        return None
    
    print(f"\n{msi_status} Therapy-specific analysis:")
    print(f"  Patients with change data: {len(ms_changes)}")
    
    return compare_treatment_changes_ms(ms_changes, msi_status, save_path)

def analyze_msi_vs_mss_programs(adata, save_path):
    """Direct comparison of MP patterns between MSI and MSS tumors"""
    
    print("\nDirect MSI vs MSS Meta-Program Comparison")
    
    # Filter for treated samples with clear microsatellite status
    treated_mask = adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx'])
    ms_clear_mask = adata.obs['Microsatellite_Status'].isin(['MSI', 'MSS'])
    
    comparison_data = adata[treated_mask & ms_clear_mask].copy()
    
    print(f"Comparison cohort: {comparison_data.n_obs:,} cells")
    print("Microsatellite distribution:")
    print(comparison_data.obs['Microsatellite_Status'].value_counts())
    
    # Overall MP distribution comparison
    overall_comparison = run_mp_comparison_analysis(
        comparison_data,
        comparison_type='msi_vs_mss_overall',
        group_col='Microsatellite_Status', 
        group1='MSS',
        group2='MSI',
        save_path=save_path,
        description="Overall MSS vs MSI Meta-Program Distribution",
        paired_analysis=False
    )
    
    # Pre-treatment specific comparison
    pre_comparison_data = comparison_data[comparison_data.obs['Treatment_Stage'] == 'Pre'].copy()
    
    if len(pre_comparison_data) > 100:
        pre_comparison = run_mp_comparison_analysis(
            pre_comparison_data,
            comparison_type='msi_vs_mss_pre',
            group_col='Microsatellite_Status',
            group1='MSS', 
            group2='MSI',
            save_path=save_path,
            description="Pre-treatment MSS vs MSI Meta-Program Patterns",
            paired_analysis=False
        )
    else:
        pre_comparison = None
    
    # Response-stratified comparison (within non-pCR and pCR separately)
    response_stratified = {}
    
    for response in ['pCR', 'non_pCR']:
        response_data = comparison_data[comparison_data.obs['Response'] == response].copy()
        
        if len(response_data) > 50:
            msi_count = (response_data.obs['Microsatellite_Status'] == 'MSI').sum()
            mss_count = (response_data.obs['Microsatellite_Status'] == 'MSS').sum()
            
            print(f"\n{response} stratified analysis: MSI={msi_count}, MSS={mss_count}")
            
            if msi_count >= 20 and mss_count >= 20:
                response_comp = run_mp_comparison_analysis(
                    response_data,
                    comparison_type=f'msi_vs_mss_{response}',
                    group_col='Microsatellite_Status',
                    group1='MSS',
                    group2='MSI', 
                    save_path=save_path,
                    description=f"MSS vs MSI in {response} patients",
                    paired_analysis=False
                )
                response_stratified[response] = response_comp
    
    return {
        'overall_comparison': overall_comparison,
        'pre_treatment_comparison': pre_comparison,
        'response_stratified': response_stratified
    }

def identify_universal_vs_specific_resistance(msi_results, mss_results, save_path):
    """Identify universal resistance mechanisms vs microsatellite-specific ones"""
    
    print("\nIdentifying Universal vs Specific Resistance Mechanisms")
    
    if not msi_results or not mss_results:
        print("⚠️ Insufficient data for universal vs specific analysis")
        return None
    
    # Extract significant MPs from each analysis
    msi_significant = extract_significant_mps(msi_results)
    mss_significant = extract_significant_mps(mss_results)
    
    # Find overlaps and differences
    universal_mps = []
    msi_specific_mps = []
    mss_specific_mps = []
    
    all_msi_mps = set(msi_significant.keys()) if msi_significant else set()
    all_mss_mps = set(mss_significant.keys()) if mss_significant else set()
    
    # Universal: significant in both MSI and MSS
    universal_mp_names = all_msi_mps & all_mss_mps
    
    # Specific: significant in only one
    msi_specific_mp_names = all_msi_mps - all_mss_mps
    mss_specific_mp_names = all_mss_mps - all_msi_mps
    
    print(f"\n🌐 Universal resistance MPs: {list(universal_mp_names)}")
    print(f"🎯 MSI-specific resistance MPs: {list(msi_specific_mp_names)}")
    print(f"🔍 MSS-specific resistance MPs: {list(mss_specific_mp_names)}")
    
    # Create summary
    universal_vs_specific_summary = {
        'universal_mps': universal_mp_names,
        'msi_specific_mps': msi_specific_mp_names,
        'mss_specific_mps': mss_specific_mp_names,
        'msi_results': msi_significant,
        'mss_results': mss_significant
    }
    
    # Save detailed results
    with open(f"{save_path}/universal_vs_specific_resistance.pkl", 'wb') as f:
        pickle.dump(universal_vs_specific_summary, f)
    
    return universal_vs_specific_summary

def extract_significant_mps(results_dict):
    """Extract significant MPs from resistance analysis results"""
    
    significant_mps = {}
    
    if not results_dict:
        return significant_mps
    
    for analysis_type, analysis_results in results_dict.items():
        if analysis_results and 'significant_mps' in analysis_results:
            for _, row in analysis_results['significant_mps'].iterrows():
                mp_name = row['mp']
                if mp_name not in significant_mps:
                    significant_mps[mp_name] = []
                
                significant_mps[mp_name].append({
                    'analysis_type': analysis_type,
                    'effect_size': row['cohens_d'],
                    'p_corrected': row['p_corrected']
                })
    
    return significant_mps

def generate_clinical_insights(msi_results, mss_results, direct_comparison, save_path):
    """Generate clinical insights from MSI vs MSS analysis"""
    
    print("\nGenerating Clinical Insights")
    
    insights = []
    
    # Response rate insights
    if msi_results and mss_results:
        insights.append("📊 Response Rate Patterns:")
        insights.append("   - Compare intrinsic resistance programs")
        insights.append("   - Identify why MSS typically doesn't respond")
        
    # Biomarker insights  
    if direct_comparison and 'overall_comparison' in direct_comparison:
        insights.append("🎯 Biomarker Implications:")
        insights.append("   - MSI-specific programs for current immunotherapy")
        insights.append("   - MSS-specific programs for future therapeutic targeting")
        
    # Therapeutic insights
    insights.append("💊 Therapeutic Insights:")
    insights.append("   - Universal targets: work across microsatellite status")
    insights.append("   - Specific targets: tailored to MSI or MSS tumors")
    
    # Save insights
    with open(f"{save_path}/clinical_insights.txt", 'w') as f:
        f.write("\n".join(insights))
    
    print("\n".join(insights))
    
    return insights

# Helper functions for MSS-specific analysis
def calculate_treatment_induced_changes_ms(ms_data, msi_status):
    """Calculate treatment changes for specific microsatellite status"""
    return calculate_treatment_induced_changes(ms_data, slice(None))  # Use existing function

def compare_treatment_changes_ms(changes_df, msi_status, save_path):
    """Compare treatment changes for specific microsatellite status"""
    # Modify the existing function to work with MS-specific data
    return compare_treatment_changes(changes_df, save_path)

# MSI-H CRC Specific Resistance Analysis Functions
def calculate_mp_fractions_per_sample(adata, mp_column='MP_assignment'):
    """
    Calculate MP fractions for each sample
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with cells
    mp_column : str
        Column name containing MP assignments
    
    Returns:
    --------
    pd.DataFrame
        Sample-level MP fractions with clinical metadata
    """
    
    print("📊 Calculating MP fractions per sample...")
    
    sample_mp_data = []
    
    for sample_id in adata.obs['Sample_ID'].unique():
        sample_cells = adata[adata.obs['Sample_ID'] == sample_id]
        
        if len(sample_cells) < 10:  # Minimum cells per sample
            continue
        
        # Get sample metadata (taking first cell's metadata as sample-level)
        sample_metadata = sample_cells.obs.iloc[0]
        
        # Calculate MP fractions
        mp_counts = sample_cells.obs[mp_column].value_counts()
        total_cells = len(sample_cells)
        
        # Get all unique MPs in the dataset
        all_mps = adata.obs[mp_column].unique()
        
        for mp in all_mps:
            if mp != 'Unresolved' and pd.notna(mp):
                fraction = mp_counts.get(mp, 0) / total_cells
                
                sample_mp_data.append({
                    'Sample_ID': sample_id,
                    'Patient_ID': sample_metadata.get('Patient_ID', sample_id),
                    'MP': mp,
                    'MP_fraction': fraction,
                    'total_cells': total_cells,
                    'mp_cell_count': mp_counts.get(mp, 0),
                    
                    # Clinical metadata
                    'Treatment_Strategy': sample_metadata.get('Treatment_Strategy', 'Unknown'),
                    'Microsatellite_Status': sample_metadata.get('Microsatellite_Status', 'Unknown'),
                    'Response': sample_metadata.get('Response', 'Unknown'),
                    'Treatment_Stage': sample_metadata.get('Treatment_Stage', 'Unknown'),
                    'Gender': sample_metadata.get('Gender', 'Unknown'),
                    'Age': sample_metadata.get('Age', np.nan),
                    'study': sample_metadata.get('study', 'Unknown'),
                    'Tissue_Type': sample_metadata.get('Tissue_Type', 'Unknown')
                })
    
    sample_df = pd.DataFrame(sample_mp_data)
    
    print(f"✅ Calculated MP fractions for {len(sample_df['Sample_ID'].unique())} samples")
    print(f"   Total MP-sample combinations: {len(sample_df)}")
    
    return sample_df

def task1_intrinsic_resistance_msi_pretreatment(adata, save_path, mp_column='MP_assignment'):
    """
    Task 1: Compare MP fractions between pCR and non-pCR in MSI pre-treatment samples
    
    Focus: Intrinsic resistance mechanisms in MSI tumors
    """
    
    print("\n🎯 Task 1: Intrinsic Resistance Analysis (MSI Pre-treatment)")
    print("="*70)
    
    # Filter criteria
    filter_mask = (
        (adata.obs['Microsatellite_Status'] == 'MSI') &
        (adata.obs['Treatment_Stage'] == 'Pre') &
        (adata.obs['Response'].isin(['pCR', 'non_pCR'])) &
        (adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx']))
    )
    
    filtered_data = adata[filter_mask].copy()
    
    print(f"📊 Data Overview:")
    print(f"   Total MSI pre-treatment cells: {filtered_data.n_obs:,}")
    print(f"   Samples: {filtered_data.obs['Sample_ID'].nunique()}")
    print(f"   Patients: {filtered_data.obs['Patient_ID'].nunique()}")
    
    # Response distribution
    response_dist = filtered_data.obs['Response'].value_counts()
    print(f"\n   Response distribution:")
    for response, count in response_dist.items():
        percentage = count / filtered_data.n_obs * 100
        print(f"     {response}: {count:,} cells ({percentage:.1f}%)")
    
    # Treatment strategy distribution
    treatment_dist = filtered_data.obs['Treatment_Strategy'].value_counts()
    print(f"\n   Treatment strategies:")
    for treatment, count in treatment_dist.items():
        print(f"     {treatment}: {count:,} cells")
    
    # Calculate sample-level MP fractions
    sample_mp_df = calculate_mp_fractions_per_sample(filtered_data, mp_column)
    
    # Statistical comparison: pCR vs non_pCR
    print(f"\n🔬 Statistical Analysis: pCR vs non_pCR")
    
    results = []
    
    for mp in sample_mp_df['MP'].unique():
        mp_data = sample_mp_df[sample_mp_df['MP'] == mp]
        
        pcr_fractions = mp_data[mp_data['Response'] == 'pCR']['MP_fraction']
        nonpcr_fractions = mp_data[mp_data['Response'] == 'non_pCR']['MP_fraction']
        
        if len(pcr_fractions) >= 3 and len(nonpcr_fractions) >= 3:
            # Perform t-test
            t_stat, p_value = stats.ttest_ind(nonpcr_fractions, pcr_fractions)
            
            # Calculate effect size (Cohen's d)
            pooled_std = np.sqrt(((len(nonpcr_fractions) - 1) * np.var(nonpcr_fractions, ddof=1) + 
                                (len(pcr_fractions) - 1) * np.var(pcr_fractions, ddof=1)) / 
                               (len(nonpcr_fractions) + len(pcr_fractions) - 2))
            
            cohens_d = (np.mean(nonpcr_fractions) - np.mean(pcr_fractions)) / pooled_std if pooled_std > 0 else 0
            
            results.append({
                'MP': mp,
                'nonpcr_mean': np.mean(nonpcr_fractions),
                'nonpcr_std': np.std(nonpcr_fractions),
                'nonpcr_n_samples': len(nonpcr_fractions),
                'pcr_mean': np.mean(pcr_fractions),
                'pcr_std': np.std(pcr_fractions),
                'pcr_n_samples': len(pcr_fractions),
                'mean_difference': np.mean(nonpcr_fractions) - np.mean(pcr_fractions),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d,
                'analysis_type': 'intrinsic_resistance'
            })
            
            print(f"\n{mp}:")
            print(f"  non_pCR: {np.mean(nonpcr_fractions):.4f} ± {np.std(nonpcr_fractions):.4f} (n={len(nonpcr_fractions)})")
            print(f"  pCR:     {np.mean(pcr_fractions):.4f} ± {np.std(pcr_fractions):.4f} (n={len(pcr_fractions)})")
            print(f"  Difference: {np.mean(nonpcr_fractions) - np.mean(pcr_fractions):.4f}")
            print(f"  t-test: t={t_stat:.3f}, p={p_value:.3e}, Cohen's d={cohens_d:.3f}")
    
    # Convert to DataFrame and apply multiple testing correction
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # FDR correction
        _, p_corrected, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        results_df['p_corrected'] = p_corrected
        
        # Identify significant intrinsic resistance MPs
        significant_mps = results_df[
            (results_df['p_corrected'] < 0.05) & 
            (abs(results_df['cohens_d']) > 0.5)
        ].copy()
        
        print(f"\n🌟 Significant Intrinsic Resistance MPs (FDR < 0.05, |Cohen's d| > 0.5):")
        if len(significant_mps) > 0:
            for _, row in significant_mps.iterrows():
                direction = "higher" if row['cohens_d'] > 0 else "lower"
                print(f"  {row['MP']}: {direction} in non_pCR")
                print(f"    Mean difference: {row['mean_difference']:.4f}")
                print(f"    Effect size: {row['cohens_d']:.3f}")
                print(f"    Corrected p-value: {row['p_corrected']:.3e}")
        else:
            print("  None found with current thresholds")
        
    # Save sample-level data for R analysis
    sample_mp_df.to_csv(f"{save_path}/task1_intrinsic_resistance_sample_data.csv", index=False)
    
    # Save statistical results
    if len(results_df) > 0:
        results_df.to_csv(f"{save_path}/task1_intrinsic_resistance_statistical_results.csv", index=False)
    
    # Save filtered cell-level data for additional analysis
    # filtered_data.obs.to_csv(f"{save_path}/task1_intrinsic_resistance_cell_metadata.csv")
    
    print(f"\n✅ Task 1 completed!")
    print(f"📁 Results saved to {save_path}/")
    print(f"   - task1_intrinsic_resistance_sample_data.csv (for R plotting)")
    print(f"   - task1_intrinsic_resistance_statistical_results.csv")
    print(f"   - task1_intrinsic_resistance_cell_metadata.csv")
    
    return {
        'sample_data': sample_mp_df,
        'statistical_results': results_df,
        'significant_mps': significant_mps if len(results_df) > 0 and len(significant_mps) > 0 else None,
        # 'filtered_adata': filtered_data
    }

def task2_acquired_resistance_msi_nonpcr(adata, save_path, mp_column='MP_assignment'):
    """
    Task 2: Compare MP fractions between Pre and Post in MSI non_pCR samples
    
    Focus: Acquired resistance mechanisms during treatment
    """
    
    print("\n🔥 Task 2: Acquired Resistance Analysis (MSI non_pCR Pre vs Post)")
    print("="*70)
    
    # Filter criteria
    filter_mask = (
        (adata.obs['Microsatellite_Status'] == 'MSI') &
        (adata.obs['Response'] == 'non_pCR') &
        (adata.obs['Treatment_Stage'].isin(['Pre', 'Post'])) &
        (adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx']))
    )
    
    filtered_data = adata[filter_mask].copy()
    
    print(f"📊 Data Overview:")
    print(f"   Total MSI non_pCR cells: {filtered_data.n_obs:,}")
    print(f"   Samples: {filtered_data.obs['Sample_ID'].nunique()}")
    print(f"   Patients: {filtered_data.obs['Patient_ID'].nunique()}")
    
    # Treatment stage distribution
    stage_dist = filtered_data.obs['Treatment_Stage'].value_counts()
    print(f"\n   Treatment stage distribution:")
    for stage, count in stage_dist.items():
        print(f"     {stage}: {count:,} cells")
    
    # Check for paired samples (same patient with both Pre and Post)
    patient_stages = filtered_data.obs.groupby('Patient_ID')['Treatment_Stage'].unique()
    paired_patients = [pid for pid, stages in patient_stages.items() 
                      if 'Pre' in stages and 'Post' in stages]
    
    print(f"\n   Patients with paired Pre/Post samples: {len(paired_patients)}")
    
    if len(paired_patients) < 3:
        print("⚠️ Warning: Very few paired samples for robust analysis")
    
    # Calculate sample-level MP fractions
    sample_mp_df = calculate_mp_fractions_per_sample(filtered_data, mp_column)
    
    # Statistical comparison: Post vs Pre
    print(f"\n🔬 Statistical Analysis: Post vs Pre")
    
    results = []
    
    for mp in sample_mp_df['MP'].unique():
        mp_data = sample_mp_df[sample_mp_df['MP'] == mp]
        
        pre_fractions = mp_data[mp_data['Treatment_Stage'] == 'Pre']['MP_fraction']
        post_fractions = mp_data[mp_data['Treatment_Stage'] == 'Post']['MP_fraction']
        
        if len(pre_fractions) >= 3 and len(post_fractions) >= 3:
            # Check if we can do paired analysis
            pre_patients = set(mp_data[mp_data['Treatment_Stage'] == 'Pre']['Patient_ID'])
            post_patients = set(mp_data[mp_data['Treatment_Stage'] == 'Post']['Patient_ID'])
            common_patients = pre_patients & post_patients
            
            if len(common_patients) >= 3:
                # Paired t-test
                paired_pre = []
                paired_post = []
                
                for patient in common_patients:
                    pre_val = mp_data[(mp_data['Patient_ID'] == patient) & 
                                    (mp_data['Treatment_Stage'] == 'Pre')]['MP_fraction']
                    post_val = mp_data[(mp_data['Patient_ID'] == patient) & 
                                     (mp_data['Treatment_Stage'] == 'Post')]['MP_fraction']
                    
                    if len(pre_val) > 0 and len(post_val) > 0:
                        paired_pre.append(pre_val.iloc[0])
                        paired_post.append(post_val.iloc[0])
                
                if len(paired_pre) >= 3:
                    t_stat, p_value = stats.ttest_rel(paired_post, paired_pre)
                    
                    # Effect size for paired data
                    differences = np.array(paired_post) - np.array(paired_pre)
                    cohens_d = np.mean(differences) / np.std(differences, ddof=1) if np.std(differences, ddof=1) > 0 else 0
                    
                    test_type = 'paired_t_test'
                    n_samples = len(paired_pre)
                    
                else:
                    continue
            else:
                # Unpaired t-test
                t_stat, p_value = stats.ttest_ind(post_fractions, pre_fractions)
                
                # Effect size
                pooled_std = np.sqrt(((len(post_fractions) - 1) * np.var(post_fractions, ddof=1) + 
                                    (len(pre_fractions) - 1) * np.var(pre_fractions, ddof=1)) / 
                                   (len(post_fractions) + len(pre_fractions) - 2))
                
                cohens_d = (np.mean(post_fractions) - np.mean(pre_fractions)) / pooled_std if pooled_std > 0 else 0
                
                test_type = 'unpaired_t_test'
                n_samples = f"Pre: {len(pre_fractions)}, Post: {len(post_fractions)}"
            
            results.append({
                'MP': mp,
                'post_mean': np.mean(post_fractions),
                'post_std': np.std(post_fractions),
                'pre_mean': np.mean(pre_fractions),
                'pre_std': np.std(pre_fractions),
                'mean_difference': np.mean(post_fractions) - np.mean(pre_fractions),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d,
                'test_type': test_type,
                'n_samples': n_samples,
                'analysis_type': 'acquired_resistance'
            })
            
            print(f"\n{mp} ({test_type}):")
            print(f"  Post: {np.mean(post_fractions):.4f} ± {np.std(post_fractions):.4f}")
            print(f"  Pre:  {np.mean(pre_fractions):.4f} ± {np.std(pre_fractions):.4f}")
            print(f"  Change: {np.mean(post_fractions) - np.mean(pre_fractions):.4f}")
            print(f"  t-test: t={t_stat:.3f}, p={p_value:.3e}, Cohen's d={cohens_d:.3f}")
    
    # Convert to DataFrame and apply multiple testing correction
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # FDR correction
        _, p_corrected, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        results_df['p_corrected'] = p_corrected
        
        # Identify significant acquired resistance MPs (increased post-treatment)
        significant_mps = results_df[
            (results_df['p_corrected'] < 0.05) & 
            (results_df['cohens_d'] > 0.5) &  # Post > Pre
            (results_df['mean_difference'] > 0)  # Actually increased
        ].copy()
        
        print(f"\n🌟 Significant Acquired Resistance MPs (FDR < 0.05, Cohen's d > 0.5, increased):")
        if len(significant_mps) > 0:
            for _, row in significant_mps.iterrows():
                print(f"  {row['MP']}: increased post-treatment")
                print(f"    Mean change: +{row['mean_difference']:.4f}")
                print(f"    Effect size: {row['cohens_d']:.3f}")
                print(f"    Corrected p-value: {row['p_corrected']:.3e}")
        else:
            print("  None found with current thresholds")
        
    # Save sample-level data for R analysis
    sample_mp_df.to_csv(f"{save_path}/task2_acquired_resistance_sample_data.csv", index=False)
    
    # Save statistical results
    if len(results_df) > 0:
        results_df.to_csv(f"{save_path}/task2_acquired_resistance_statistical_results.csv", index=False)
    
    # Save paired patients info
    paired_info = pd.DataFrame({
        'Patient_ID': paired_patients,
        'has_paired_samples': True
    })
    paired_info.to_csv(f"{save_path}/task2_paired_patients_info.csv", index=False)
    
    # Save filtered cell-level data
    filtered_data.obs.to_csv(f"{save_path}/task2_acquired_resistance_cell_metadata.csv")
    
    print(f"\n✅ Task 2 completed!")
    print(f"📁 Results saved to {save_path}/")
    print(f"   - task2_acquired_resistance_sample_data.csv (for R plotting)")
    print(f"   - task2_acquired_resistance_statistical_results.csv")
    print(f"   - task2_paired_patients_info.csv")
    print(f"   - task2_acquired_resistance_cell_metadata.csv")
    
    return {
        'sample_data': sample_mp_df,
        'statistical_results': results_df,
        'significant_mps': significant_mps if len(results_df) > 0 and len(significant_mps) > 0 else None,
        'paired_patients': paired_patients,
        'filtered_adata': filtered_data
    }

def task3_msi_mss_mp_similarity(adata, save_path, mp_column='MP_assignment'):
    """
    Task 3: Analyze MP similarity between MSI and MSS tumors
    
    Focus: Understanding shared vs distinct MP patterns
    """
    
    print("\n⚡ Task 3: MSI vs MSS Meta-Program Similarity Analysis")
    print("="*70)
    
    # Filter for treated samples with clear microsatellite status
    filter_mask = (
        (adata.obs['Microsatellite_Status'].isin(['MSI', 'MSS'])) &
        (adata.obs['Treatment_Stage'].isin(['Pre'])) &
        (adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx']))
    )
    
    filtered_data = adata[filter_mask].copy()
    
    print(f"📊 Data Overview:")
    print(f"   Total treated cells: {filtered_data.n_obs:,}")
    print(f"   Samples: {filtered_data.obs['Sample_ID'].nunique()}")
    print(f"   Patients: {filtered_data.obs['Patient_ID'].nunique()}")
    
    # Microsatellite status distribution
    ms_dist = filtered_data.obs['Microsatellite_Status'].value_counts()
    print(f"\n   Microsatellite status distribution:")
    for status, count in ms_dist.items():
        percentage = count / filtered_data.n_obs * 100
        print(f"     {status}: {count:,} cells ({percentage:.1f}%)")
    
    # Calculate sample-level MP fractions
    sample_mp_df = calculate_mp_fractions_per_sample(filtered_data, mp_column)
    
    # Statistical comparison: MSS vs MSI
    print(f"\n🔬 Statistical Analysis: MSS vs MSI")
    
    results = []
    
    for mp in sample_mp_df['MP'].unique():
        mp_data = sample_mp_df[sample_mp_df['MP'] == mp]
        
        msi_fractions = mp_data[mp_data['Microsatellite_Status'] == 'MSI']['MP_fraction']
        mss_fractions = mp_data[mp_data['Microsatellite_Status'] == 'MSS']['MP_fraction']
        
        if len(msi_fractions) >= 3 and len(mss_fractions) >= 3:
            # Perform t-test
            t_stat, p_value = stats.ttest_ind(mss_fractions, msi_fractions)
            
            # Calculate effect size (Cohen's d)
            pooled_std = np.sqrt(((len(mss_fractions) - 1) * np.var(mss_fractions, ddof=1) + 
                                (len(msi_fractions) - 1) * np.var(msi_fractions, ddof=1)) / 
                               (len(mss_fractions) + len(msi_fractions) - 2))
            
            cohens_d = (np.mean(mss_fractions) - np.mean(msi_fractions)) / pooled_std if pooled_std > 0 else 0
            
            results.append({
                'MP': mp,
                'mss_mean': np.mean(mss_fractions),
                'mss_std': np.std(mss_fractions),
                'mss_n_samples': len(mss_fractions),
                'msi_mean': np.mean(msi_fractions),
                'msi_std': np.std(msi_fractions),
                'msi_n_samples': len(msi_fractions),
                'mean_difference': np.mean(mss_fractions) - np.mean(msi_fractions),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d,
                'analysis_type': 'msi_mss_comparison'
            })
            
            print(f"\n{mp}:")
            print(f"  MSS: {np.mean(mss_fractions):.4f} ± {np.std(mss_fractions):.4f} (n={len(mss_fractions)})")
            print(f"  MSI: {np.mean(msi_fractions):.4f} ± {np.std(msi_fractions):.4f} (n={len(msi_fractions)})")
            print(f"  Difference: {np.mean(mss_fractions) - np.mean(msi_fractions):.4f}")
            print(f"  t-test: t={t_stat:.3f}, p={p_value:.3e}, Cohen's d={cohens_d:.3f}")
    
    # Convert to DataFrame and apply multiple testing correction
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # FDR correction
        _, p_corrected, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        results_df['p_corrected'] = p_corrected
        
        # Classify MPs based on differential patterns
        significant_mps = results_df[
            (results_df['p_corrected'] < 0.05) & 
            (abs(results_df['cohens_d']) > 0.5)
        ].copy()
        
        # Categorize MPs
        mss_enriched = significant_mps[significant_mps['cohens_d'] > 0.5]  # MSS > MSI
        msi_enriched = significant_mps[significant_mps['cohens_d'] < -0.5]  # MSI > MSS
        shared_mps = results_df[
            (results_df['p_corrected'] >= 0.05) |  # Not significantly different
            (abs(results_df['cohens_d']) <= 0.2)   # Small effect size
        ]
        
        print(f"\n🌟 MP Classification Results:")
        print(f"  MSS-enriched MPs: {len(mss_enriched)}")
        if len(mss_enriched) > 0:
            for _, row in mss_enriched.iterrows():
                print(f"    {row['MP']}: {row['cohens_d']:.3f} (p={row['p_corrected']:.3e})")
        
        print(f"  MSI-enriched MPs: {len(msi_enriched)}")
        if len(msi_enriched) > 0:
            for _, row in msi_enriched.iterrows():
                print(f"    {row['MP']}: {row['cohens_d']:.3f} (p={row['p_corrected']:.3e})")
        
        print(f"  Shared MPs (similar between MSI/MSS): {len(shared_mps)}")
        if len(shared_mps) > 0:
            for _, row in shared_mps.head(5).iterrows():  # Show top 5
                print(f"    {row['MP']}: {row['cohens_d']:.3f} (p={row['p_corrected']:.3e})")
    
    # Calculate MP correlation matrix between MSI and MSS
    print(f"\n📊 Calculating MP correlation between MSI and MSS...")
    
    # Create pivot table: samples x MPs
    pivot_df = sample_mp_df.pivot_table(
        index='Sample_ID', 
        columns='MP', 
        values='MP_fraction', 
        fill_value=0
    )
    
    # Add microsatellite status
    sample_ms_status = sample_mp_df.groupby('Sample_ID')['Microsatellite_Status'].first()
    pivot_df['Microsatellite_Status'] = sample_ms_status
    
    # Calculate mean MP fractions for MSI and MSS
    msi_means = pivot_df[pivot_df['Microsatellite_Status'] == 'MSI'].drop('Microsatellite_Status', axis=1).mean()
    mss_means = pivot_df[pivot_df['Microsatellite_Status'] == 'MSS'].drop('Microsatellite_Status', axis=1).mean()
    
    # Calculate correlation
    mp_correlation = stats.pearsonr(msi_means, mss_means)
    
    print(f"   MP pattern correlation (MSI vs MSS): r={mp_correlation[0]:.3f}, p={mp_correlation[1]:.3e}")
    
    # Create similarity summary
    similarity_summary = {
        'mp_correlation': mp_correlation[0],
        'mp_correlation_pvalue': mp_correlation[1],
        'total_mps': len(results_df),
        'mss_enriched_mps': len(mss_enriched) if len(results_df) > 0 else 0,
        'msi_enriched_mps': len(msi_enriched) if len(results_df) > 0 else 0,
        'shared_mps': len(shared_mps) if len(results_df) > 0 else 0,
        'msi_mean_fractions': msi_means.to_dict(),
        'mss_mean_fractions': mss_means.to_dict()
    }
    
    # Save sample-level data for R analysis
    sample_mp_df.to_csv(f"{save_path}/task3_msi_mss_similarity_sample_data.csv", index=False)
    
    # Save statistical results
    if len(results_df) > 0:
        results_df.to_csv(f"{save_path}/task3_msi_mss_similarity_statistical_results.csv", index=False)
    
        # Save MP classification
        mp_classification = pd.DataFrame({
            'MP': results_df['MP'],
            'category': pd.cut(results_df['cohens_d'], 
                             bins=[-np.inf, -0.5, -0.2, 0.2, 0.5, np.inf],
                             labels=['MSI_enriched', 'MSI_slight', 'Shared', 'MSS_slight', 'MSS_enriched']),
            'cohens_d': results_df['cohens_d'],
            'p_corrected': results_df['p_corrected'],
            'significant': results_df['p_corrected'] < 0.05,
            'mss_mean': results_df['mss_mean'],
            'msi_mean': results_df['msi_mean'],
            'mean_difference': results_df['mean_difference']
        })
        mp_classification.to_csv(f"{save_path}/task3_mp_classification.csv", index=False)
    
    # Save pivot table for correlation analysis
    pivot_df.to_csv(f"{save_path}/task3_mp_pivot_table.csv")
    
    # Save similarity summary
    with open(f"{save_path}/task3_similarity_summary.json", 'w') as f:
        json.dump(similarity_summary, f, indent=2)
    
    # Save filtered cell-level data
    filtered_data.obs.to_csv(f"{save_path}/task3_msi_mss_cell_metadata.csv")
    
    print(f"\n✅ Task 3 completed!")
    print(f"📁 Results saved to {save_path}/")
    print(f"   - task3_msi_mss_similarity_sample_data.csv (for R plotting)")
    print(f"   - task3_msi_mss_similarity_statistical_results.csv")
    print(f"   - task3_mp_classification.csv")
    print(f"   - task3_mp_pivot_table.csv")
    print(f"   - task3_similarity_summary.json")
    print(f"   - task3_msi_mss_cell_metadata.csv")
    
    return {
        'sample_data': sample_mp_df,
        'statistical_results': results_df,
        'mp_classification': mp_classification if len(results_df) > 0 else None,
        'similarity_summary': similarity_summary,
        'pivot_table': pivot_df,
        'mss_enriched_mps': mss_enriched if len(results_df) > 0 else None,
        'msi_enriched_mps': msi_enriched if len(results_df) > 0 else None,
        'shared_mps': shared_mps if len(results_df) > 0 else None,
        'filtered_adata': filtered_data
    }

# Missing Key Figures Implementation
def create_nmf_robustness_plots(robust_results, save_path):
    """Figure 1: NMF Parameter Sensitivity Analysis"""
    
    print("📊 Creating NMF robustness plots...")
    
    # Extract robustness statistics
    all_programs = robust_results['all_programs']
    robust_programs = robust_results['robust_programs']
    parameters = robust_results['parameters']
    
    # Analyze program survival by K value
    k_survival = {}
    for k in parameters['K_range']:
        total_programs = len([p for p in all_programs if p['k_value'] == k])
        surviving_programs = len([p for p in robust_programs if p['k_value'] == k])
        k_survival[k] = {'total': total_programs, 'surviving': surviving_programs, 
                        'survival_rate': surviving_programs / total_programs if total_programs > 0 else 0}
    
    # Create subplot figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Program survival by K
    k_values = list(k_survival.keys())
    survival_rates = [k_survival[k]['survival_rate'] for k in k_values]
    total_programs = [k_survival[k]['total'] for k in k_values]
    
    axes[0,0].bar(k_values, survival_rates, alpha=0.7, color='steelblue')
    axes[0,0].set_xlabel('Number of NMF Components (K)')
    axes[0,0].set_ylabel('Program Survival Rate')
    axes[0,0].set_title('NMF Robustness: Program Survival by K Value')
    axes[0,0].set_ylim(0, 1)
    
    # Add text annotations
    for i, (k, rate) in enumerate(zip(k_values, survival_rates)):
        axes[0,0].text(k, rate + 0.02, f'{rate:.2f}', ha='center', fontweight='bold')
    
    # Plot 2: Total programs generated by K
    axes[0,1].bar(k_values, total_programs, alpha=0.7, color='orange')
    axes[0,1].set_xlabel('Number of NMF Components (K)')
    axes[0,1].set_ylabel('Total Programs Generated')
    axes[0,1].set_title('Program Generation by K Value')
    
    # Plot 3: Patient representation in robust programs
    patient_counts = {}
    for program in robust_programs:
        patient = program['patient_id']
        patient_counts[patient] = patient_counts.get(patient, 0) + 1
    
    patient_program_counts = list(patient_counts.values())
    axes[1,0].hist(patient_program_counts, bins=max(10, len(set(patient_program_counts))), 
                   alpha=0.7, color='green', edgecolor='black')
    axes[1,0].set_xlabel('Number of Robust Programs per Patient')
    axes[1,0].set_ylabel('Number of Patients')
    axes[1,0].set_title('Patient Contribution to Robust Programs')
    
    # Plot 4: Robustness criteria statistics
    criteria_names = ['Total Programs', 'Cross-K Stable', 'Cross-Patient', 'Final Robust']
    criteria_counts = [
        len(all_programs),
        len(all_programs),  # This would need to be tracked from the actual function
        len(all_programs),  # This would need to be tracked from the actual function  
        len(robust_programs)
    ]
    
    axes[1,1].bar(range(len(criteria_names)), criteria_counts, alpha=0.7, color='purple')
    axes[1,1].set_xticks(range(len(criteria_names)))
    axes[1,1].set_xticklabels(criteria_names, rotation=45, ha='right')
    axes[1,1].set_ylabel('Number of Programs')
    axes[1,1].set_title('Program Filtering Pipeline')
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/missing_fig1_nmf_robustness.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ NMF robustness plots saved")

def create_mp_assignment_quality_plots(malignant_cells, save_path):
    """Figure 2: MP Assignment Quality Assessment"""
    
    print("📊 Creating MP assignment quality plots...")
    
    # Get MP score columns
    mp_score_cols = [col for col in malignant_cells.obs.columns if 'MP' in col and 'score' in col and 'normalized' not in col]
    
    if len(mp_score_cols) == 0:
        print("⚠️ No MP scores found")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: MP score distributions
    score_data = malignant_cells.obs[mp_score_cols].values
    
    axes[0,0].boxplot(score_data, labels=[col.replace('_score', '') for col in mp_score_cols])
    axes[0,0].set_xlabel('Meta-Programs')
    axes[0,0].set_ylabel('MP Score')
    axes[0,0].set_title('MP Score Distributions')
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Assignment confidence (max vs second-max score)
    max_scores = np.max(score_data, axis=1)
    second_max_scores = np.partition(score_data, -2, axis=1)[:, -2]
    confidence_scores = max_scores - second_max_scores
    
    axes[0,1].hist(confidence_scores, bins=50, alpha=0.7, color='orange', edgecolor='black')
    axes[0,1].axvline(x=np.median(confidence_scores), color='red', linestyle='--', 
                      label=f'Median: {np.median(confidence_scores):.3f}')
    axes[0,1].set_xlabel('Assignment Confidence (Max - Second Max Score)')
    axes[0,1].set_ylabel('Number of Cells')
    axes[0,1].set_title('MP Assignment Confidence Distribution')
    axes[0,1].legend()
    
    # Plot 3: MP assignment proportions
    if 'MP_assignment' in malignant_cells.obs.columns:
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        
        # Create pie chart
        axes[1,0].pie(assignment_counts.values, labels=assignment_counts.index, autopct='%1.1f%%', startangle=90)
        axes[1,0].set_title('MP Assignment Distribution')
    
    # Plot 4: Score correlation matrix
    score_corr = malignant_cells.obs[mp_score_cols].corr()
    
    im = axes[1,1].imshow(score_corr, cmap='RdBu_r', vmin=-1, vmax=1)
    axes[1,1].set_xticks(range(len(mp_score_cols)))
    axes[1,1].set_yticks(range(len(mp_score_cols)))
    axes[1,1].set_xticklabels([col.replace('_score', '') for col in mp_score_cols], rotation=45)
    axes[1,1].set_yticklabels([col.replace('_score', '') for col in mp_score_cols])
    axes[1,1].set_title('MP Score Correlation Matrix')
    
    # Add correlation values as text
    for i in range(len(mp_score_cols)):
        for j in range(len(mp_score_cols)):
            text = axes[1,1].text(j, i, f'{score_corr.iloc[i, j]:.2f}', 
                                 ha="center", va="center", color="white" if abs(score_corr.iloc[i, j]) > 0.5 else "black")
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=axes[1,1], shrink=0.8)
    cbar.set_label('Correlation Coefficient')
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/missing_fig2_mp_assignment_quality.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ MP assignment quality plots saved")

def create_clustering_dendrogram(robust_results, save_path):
    """Figure 3: Hierarchical Clustering Dendrogram"""
    
    print("📊 Creating clustering dendrogram...")
    
    if 'clustering_info' not in robust_results:
        print("⚠️ No clustering info found")
        return
    
    clustering_info = robust_results['clustering_info']
    linkage_matrix = clustering_info['linkage_matrix']
    cluster_assignments = clustering_info['cluster_assignments']
    
    # Create dendrogram
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))
    
    # Plot 1: Full dendrogram
    dendrogram(linkage_matrix, ax=ax1, color_threshold=0.7, above_threshold_color='gray')
    ax1.set_title('Hierarchical Clustering of Robust Programs')
    ax1.set_xlabel('Program Index')
    ax1.set_ylabel('Distance')
    
    # Plot 2: Cluster size distribution
    unique_clusters, cluster_sizes = np.unique(cluster_assignments, return_counts=True)
    
    bars = ax2.bar(range(len(unique_clusters)), cluster_sizes, alpha=0.7, color='steelblue')
    ax2.set_xlabel('Meta-Program ID')
    ax2.set_ylabel('Number of Programs')
    ax2.set_title('Programs per Meta-Program')
    ax2.set_xticks(range(len(unique_clusters)))
    ax2.set_xticklabels([f'MP{i}' for i in unique_clusters])
    
    # Add value labels on bars
    for bar, size in zip(bars, cluster_sizes):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1, f'{size}', 
                ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/missing_fig3_clustering_dendrogram.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ Clustering dendrogram saved")

def create_pathway_enrichment_heatmap(robust_results, save_path):
    """Figure 5: Pathway Enrichment Analysis for MPs"""
    
    print("📊 Creating pathway enrichment analysis...")
    
    if 'meta_programs' not in robust_results:
        print("⚠️ No meta-programs found")
        return
    
    meta_programs = robust_results['meta_programs']
    
    # Perform pathway enrichment for each MP
    enrichment_results = {}
    
    for mp_id, mp_data in meta_programs.items():
        signature_genes = mp_data['signature_genes']
        
        try:
            # GSEA using KEGG pathways
            enr = gp.enrichr(
                gene_list=signature_genes,
                gene_sets=['KEGG_2021_Human'],
                organism='Human',
                outdir=None,
                no_plot=True
            )
            
            # Get top 10 pathways
            top_pathways = enr.results.head(10)
            enrichment_results[f'MP{mp_id}'] = top_pathways[['Term', 'Adjusted P-value', 'Combined Score']]
            
        except Exception as e:
            print(f"⚠️ Enrichment failed for MP{mp_id}: {e}")
            continue
    
    if enrichment_results:
        # Create enrichment heatmap
        all_pathways = set()
        for mp_results in enrichment_results.values():
            all_pathways.update(mp_results['Term'].tolist())
        
        # Create matrix
        enrichment_matrix = pd.DataFrame(index=list(all_pathways), columns=list(enrichment_results.keys()))
        
        for mp, results in enrichment_results.items():
            for _, row in results.iterrows():
                pathway = row['Term']
                score = -np.log10(row['Adjusted P-value']) if row['Adjusted P-value'] > 0 else 10
                enrichment_matrix.loc[pathway, mp] = score
        
        enrichment_matrix = enrichment_matrix.fillna(0).astype(float)
        
        # Plot heatmap
        plt.figure(figsize=(12, 16))
        sns.heatmap(enrichment_matrix, cmap='Reds', cbar_kws={'label': '-log10(Adj P-value)'})
        plt.title('Pathway Enrichment Across Meta-Programs')
        plt.xlabel('Meta-Programs')
        plt.ylabel('KEGG Pathways')
        plt.tight_layout()
        plt.savefig(f"{save_path}/missing_fig5_pathway_enrichment.pdf", dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✅ Pathway enrichment heatmap saved")
    else:
        print("⚠️ No enrichment results to plot")

# Export function for R integration
def export_data_for_r_figures(malignant_cells, robust_results, save_path):
    """Export data needed for R figure generation"""
    
    print("📤 Exporting data for R figure generation...")
    
    # 1. Sample-level MP fractions for stacked barplot
    sample_mp_fractions = []
    
    for sample_id in malignant_cells.obs['Sample_ID'].unique():
        sample_cells = malignant_cells[malignant_cells.obs['Sample_ID'] == sample_id]
        sample_metadata = sample_cells.obs.iloc[0]
        
        if 'MP_assignment' in sample_cells.obs.columns:
            mp_counts = sample_cells.obs['MP_assignment'].value_counts()
            total_cells = len(sample_cells)
            
            for mp in malignant_cells.obs['MP_assignment'].unique():
                if mp != 'Unresolved':
                    fraction = mp_counts.get(mp, 0) / total_cells
                    
                    sample_mp_fractions.append({
                        'Sample_ID': sample_id,
                        'Patient_ID': sample_metadata.get('Patient_ID'),
                        'MP': mp,
                        'MP_fraction': fraction,
                        'MP_count': mp_counts.get(mp, 0),
                        'total_cells': total_cells,
                        'Microsatellite_Status': sample_metadata.get('Microsatellite_Status'),
                        'Response': sample_metadata.get('Response'),
                        'Treatment_Stage': sample_metadata.get('Treatment_Stage'),
                        'Treatment_Strategy': sample_metadata.get('Treatment_Strategy')
                    })
    
    sample_mp_df = pd.DataFrame(sample_mp_fractions)
    sample_mp_df.to_csv(f"{save_path}/r_sample_mp_fractions.csv", index=False)
    
    # 2. MP signature genes for heatmap
    if 'meta_programs' in robust_results:
        mp_signatures = {}
        all_signature_genes = set()
        
        for mp_id, mp_data in robust_results['meta_programs'].items():
            signature_genes = mp_data['signature_genes']
            mp_signatures[f'MP{mp_id}'] = signature_genes
            all_signature_genes.update(signature_genes)
        
        # Create binary signature matrix
        signature_matrix = pd.DataFrame(0, index=list(all_signature_genes), 
                                       columns=list(mp_signatures.keys()))
        
        for mp_name, genes in mp_signatures.items():
            for gene in genes:
                signature_matrix.loc[gene, mp_name] = 1
        
        signature_matrix.to_csv(f"{save_path}/r_mp_signature_matrix.csv")
        
        # Expression matrix for signature genes
        if all_signature_genes:
            available_genes = [g for g in all_signature_genes if g in malignant_cells.var_names]
            if available_genes:
                expr_subset = malignant_cells[:, available_genes]
                
                # Calculate mean expression per sample per gene
                sample_gene_expr = []
                for sample_id in malignant_cells.obs['Sample_ID'].unique():
                    sample_cells = expr_subset[expr_subset.obs['Sample_ID'] == sample_id]
                    mean_expr = sample_cells.X.mean(axis=0)
                    
                    if hasattr(mean_expr, 'A1'):  # Handle sparse matrices
                        mean_expr = mean_expr.A1
                    
                    for i, gene in enumerate(available_genes):
                        sample_gene_expr.append({
                            'Sample_ID': sample_id,
                            'Gene': gene,
                            'Mean_Expression': mean_expr[i]
                        })
                
                sample_expr_df = pd.DataFrame(sample_gene_expr)
                sample_expr_df.to_csv(f"{save_path}/r_sample_gene_expression.csv", index=False)
    
    # 3. MP similarity matrix
    if 'clustering_info' in robust_results:
        similarity_matrix = robust_results['clustering_info']['similarity_matrix']
        mp_names = [f'MP{i+1}' for i in range(len(similarity_matrix))]
        
        similarity_df = pd.DataFrame(similarity_matrix, index=mp_names, columns=mp_names)
        similarity_df.to_csv(f"{save_path}/r_mp_similarity_matrix.csv")
    
    print("✅ Data exported for R figure generation")
    print(f"📁 Files created:")
    print(f"   - r_sample_mp_fractions.csv (for stacked barplot)")
    print(f"   - r_mp_signature_matrix.csv (for signature heatmap)")
    print(f"   - r_sample_gene_expression.csv (for expression heatmap)")
    print(f"   - r_mp_similarity_matrix.csv (for MP correlation)")

    """Diagnose why clustering is problematic"""
    
    print("🔍 Diagnosing clustering issues...")
    
    robust_programs = robust_results.get('robust_programs', [])
    clustering_info = robust_results.get('clustering_info', {})
    
    if not robust_programs:
        print("❌ No robust programs found!")
        return
    
    print(f"📊 Analysis of {len(robust_programs)} robust programs:")
    
    # 1. Check program similarity distribution
    similarities = []
    for i in range(len(robust_programs)):
        for j in range(i+1, len(robust_programs)):
            genes1 = set(robust_programs[i]['genes'])
            genes2 = set(robust_programs[j]['genes'])
            jaccard = len(genes1 & genes2) / len(genes1 | genes2)
            similarities.append(jaccard)
    
    print(f"   Similarity stats:")
    print(f"   - Mean: {np.mean(similarities):.3f}")
    print(f"   - Median: {np.median(similarities):.3f}")
    print(f"   - Min: {np.min(similarities):.3f}")
    print(f"   - Max: {np.max(similarities):.3f}")
    print(f"   - 90th percentile: {np.percentile(similarities, 90):.3f}")
    
    # 2. Check if programs are too similar (causing massive clusters)
    high_similarity_count = sum(1 for s in similarities if s > 0.5)
    print(f"   - High similarity pairs (>0.5): {high_similarity_count}/{len(similarities)}")
    
    # 3. Patient distribution
    patient_counts = {}
    for program in robust_programs:
        patient = program.get('patient_id', 'Unknown')
        patient_counts[patient] = patient_counts.get(patient, 0) + 1
    
    print(f"   Programs per patient: {dict(patient_counts)}")
    
    # 4. K-factor distribution
    k_counts = {}
    for program in robust_programs:
        k = program.get('K', 'Unknown')
        k_counts[k] = k_counts.get(k, 0) + 1
    
    print(f"   Programs per K: {dict(k_counts)}")
    
    return {
        'similarities': similarities,
        'patient_counts': patient_counts,
        'k_counts': k_counts
    }

from scipy.stats import zscore

def calculate_mp_sample_distribution(adata, mp_names):
    """Calculate how many samples each MP appears in"""
    
    distribution = {}
    mp_names = ["MP"+ str(x) for x in mp_names]
    
    for mp_name in mp_names:
        if 'MP_assignment' in adata.obs.columns and 'Sample_ID' in adata.obs.columns:
            # Count samples where this MP is present
            mp_cells = adata.obs['MP_assignment'] == mp_name
            samples_with_mp = adata.obs.loc[mp_cells, 'Sample_ID'].nunique()
            distribution[mp_name] = samples_with_mp
        else:
            distribution[mp_name] = 0
    
    return distribution