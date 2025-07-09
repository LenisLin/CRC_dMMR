# =============================================================================
# R CODE: Hierarchical Clustering Heatmap (Figure 1D Style)
# Recreating the CSCC paper's hierarchical clustering of NMF programs
# =============================================================================

create_hierarchical_mp_heatmap <- function(patient_data, mp_signatures, output_dir, figures_dir) {
  """
  Create hierarchical clustering heatmap similar to Figure 1D in CSCC paper
  
  Parameters:
  -----------
  patient_data : data.frame
    Patient-level data with MP proportions
  mp_signatures : data.frame  
    MP signature genes matrix
  output_dir : str
    Output directory for results
  figures_dir : str
    Directory for saving figures
  """
  
  cat("🔥 Creating Hierarchical Clustering Heatmap (Figure 1D style)...\n")
  
  # ========================================================================
  # 1. CALCULATE JACCARD SIMILARITIES BETWEEN MP PROGRAMS
  # ========================================================================
  
  cat("  📊 Calculating Jaccard similarities between MP programs...\n")
  
  # Get MP signature genes (excluding the gene name column)
  mp_cols <- names(mp_signatures)[!names(mp_signatures) %in% c("V1", "gene", "Gene")]
  
  # Create MP signature lists
  mp_gene_lists <- list()
  for (mp_col in mp_cols) {
    # Get genes where this MP has value 1
    mp_genes <- mp_signatures %>%
      filter(.data[[mp_col]] == 1) %>%
      pull(if("gene" %in% names(.)) gene else V1)
    
    mp_gene_lists[[mp_col]] <- mp_genes
  }
  
  # Calculate Jaccard similarities between all MP pairs
  calculate_jaccard <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    if (union == 0) return(0)
    return(intersection / union)
  }
  
  # Create similarity matrix
  n_mps <- length(mp_gene_lists)
  mp_names <- names(mp_gene_lists)
  jaccard_matrix <- matrix(0, nrow = n_mps, ncol = n_mps)
  rownames(jaccard_matrix) <- mp_names
  colnames(jaccard_matrix) <- mp_names
  
  # Fill similarity matrix
  for (i in 1:n_mps) {
    for (j in 1:n_mps) {
      jaccard_matrix[i, j] <- calculate_jaccard(mp_gene_lists[[i]], mp_gene_lists[[j]])
    }
  }
  
  cat("  ✅ Jaccard similarity matrix calculated\n")
  cat("     Matrix dimensions:", dim(jaccard_matrix), "\n")
  cat("     Average similarity:", round(mean(jaccard_matrix[jaccard_matrix != 1]), 3), "\n")
  
  # ========================================================================
  # 2. PREPARE PATIENT-LEVEL MP PROPORTION DATA
  # ========================================================================
  
  cat("  👥 Preparing patient-level MP proportion data...\n")
  
  # Get MP proportion columns
  mp_prop_cols <- names(patient_data)[grepl("_proportion$", names(patient_data))]
  
  # Create patient MP matrix
  patient_mp_matrix <- patient_data %>%
    select(Patient_ID, all_of(mp_prop_cols)) %>%
    column_to_rownames("Patient_ID") %>%
    as.matrix()
  
  # Clean column names (remove _proportion suffix)
  colnames(patient_mp_matrix) <- str_remove(colnames(patient_mp_matrix), "_proportion")
  
  # Filter out patients with too many missing values
  patient_mp_matrix <- patient_mp_matrix[rowSums(is.na(patient_mp_matrix)) < ncol(patient_mp_matrix) * 0.5, ]
  
  # Replace remaining NAs with 0
  patient_mp_matrix[is.na(patient_mp_matrix)] <- 0
  
  cat("     Patient matrix dimensions:", dim(patient_mp_matrix), "\n")
  cat("     Patients included:", nrow(patient_mp_matrix), "\n")
  
  # ========================================================================
  # 3. CREATE PATIENT ANNOTATIONS
  # ========================================================================
  
  cat("  🏷️ Creating patient annotations...\n")
  
  # Prepare patient annotations
  patient_annotations <- patient_data %>%
    filter(Patient_ID %in% rownames(patient_mp_matrix)) %>%
    select(Patient_ID, Treatment_Strategy, Microsatellite_Status, Response, study) %>%
    column_to_rownames("Patient_ID")
  
  # Clean up annotation data
  patient_annotations <- patient_annotations %>%
    mutate(
      Treatment_Strategy = ifelse(is.na(Treatment_Strategy) | Treatment_Strategy == "", 
                                 "Unknown", Treatment_Strategy),
      Microsatellite_Status = ifelse(is.na(Microsatellite_Status) | Microsatellite_Status == "", 
                                   "Unknown", Microsatellite_Status),
      Response = ifelse(is.na(Response) | Response == "", "Unknown", Response),
      study = ifelse(is.na(study) | study == "", "Unknown", study)
    )
  
  # ========================================================================
  # 4. PERFORM HIERARCHICAL CLUSTERING
  # ========================================================================
  
  cat("  🌳 Performing hierarchical clustering...\n")
  
  # Hierarchical clustering for MPs (using 1 - Jaccard similarity as distance)
  mp_dist <- as.dist(1 - jaccard_matrix)
  mp_hclust <- hclust(mp_dist, method = "ward.D2")
  
  # Hierarchical clustering for patients (using correlation distance)
  patient_dist <- dist(patient_mp_matrix, method = "euclidean")
  patient_hclust <- hclust(patient_dist, method = "ward.D2")
  
  # ========================================================================
  # 5. CREATE COLOR SCHEMES
  # ========================================================================
  
  cat("  🎨 Setting up color schemes...\n")
  
  # Color scheme for treatment strategies
  treatment_colors <- c(
    "Anti-PD1" = "#E74C3C",
    "Anti-PD1 plus Celecoxib" = "#3498DB", 
    "Anti-PD1 plus CapeOx" = "#F39C12",
    "Unknown" = "#95A5A6"
  )
  
  # Color scheme for MSI status
  msi_colors <- c(
    "MSI" = "#E74C3C",
    "MSS" = "#3498DB",
    "Unknown" = "#95A5A6"
  )
  
  # Color scheme for response
  response_colors <- c(
    "pCR" = "#27AE60",
    "non_pCR" = "#E67E22", 
    "Unknown" = "#95A5A6"
  )
  
  # Color scheme for studies
  studies <- unique(patient_annotations$study)
  study_colors <- setNames(rainbow(length(studies)), studies)
  
  # Create annotation object for ComplexHeatmap
  col_annotation <- HeatmapAnnotation(
    Treatment = patient_annotations$Treatment_Strategy,
    MSI_Status = patient_annotations$Microsatellite_Status,
    Response = patient_annotations$Response,
    Study = patient_annotations$study,
    col = list(
      Treatment = treatment_colors[names(treatment_colors) %in% unique(patient_annotations$Treatment_Strategy)],
      MSI_Status = msi_colors[names(msi_colors) %in% unique(patient_annotations$Microsatellite_Status)],
      Response = response_colors[names(response_colors) %in% unique(patient_annotations$Response)],
      Study = study_colors[names(study_colors) %in% unique(patient_annotations$study)]
    ),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
      Treatment = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)),
      MSI_Status = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)),
      Response = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)),
      Study = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8))
    )
  )
  
  # ========================================================================
  # 6. CREATE THE MAIN HEATMAP (PATIENT MP PROPORTIONS)
  # ========================================================================
  
  cat("  🔥 Creating main heatmap...\n")
  
  # Create the main heatmap with patient MP proportions
  main_heatmap <- Heatmap(
    t(patient_mp_matrix),  # Transpose so MPs are rows, patients are columns
    name = "MP\nProportion",
    
    # Clustering
    cluster_rows = mp_hclust,
    cluster_columns = patient_hclust,
    
    # Colors
    col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    
    # Annotations
    top_annotation = col_annotation,
    
    # Appearance
    show_row_names = TRUE,
    show_column_names = FALSE,  # Too many patients to show names
    row_names_gp = gpar(fontsize = 10),
    
    # Heatmap body
    rect_gp = gpar(col = "white", lwd = 0.5),
    
    # Title
    column_title = "Patients",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title = "Meta-Programs", 
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    # Legend
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8),
      legend_direction = "vertical",
      legend_width = unit(4, "cm")
    )
  )
  
  # ========================================================================
  # 7. CREATE JACCARD SIMILARITY HEATMAP (SIDE PANEL)
  # ========================================================================
  
  cat("  📊 Creating Jaccard similarity heatmap...\n")
  
  # Create Jaccard similarity heatmap
  jaccard_heatmap <- Heatmap(
    jaccard_matrix,
    name = "Jaccard\nSimilarity",
    
    # Clustering (same as main heatmap for rows)
    cluster_rows = mp_hclust,
    cluster_columns = mp_hclust,
    
    # Colors
    col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),
    
    # Appearance
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    
    # Cell annotations (show values)
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", jaccard_matrix[i, j]), x, y, 
               gp = gpar(fontsize = 6))
    },
    
    # Title
    column_title = "MP Similarity",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    
    # Size
    width = unit(4, "cm"),
    height = unit(4, "cm"),
    
    # Legend
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 6),
      legend_direction = "vertical"
    )
  )
  
  # ========================================================================
  # 8. COMBINE HEATMAPS AND SAVE
  # ========================================================================
  
  cat("  💾 Combining heatmaps and saving...\n")
  
  # Combine heatmaps
  combined_heatmap <- main_heatmap + jaccard_heatmap
  
  # Save high-quality PDF
  pdf(file.path(figures_dir, "hierarchical_mp_clustering_heatmap.pdf"), 
      width = 16, height = 10)
  draw(combined_heatmap, 
       column_title = "Meta-Program Hierarchical Clustering Analysis",
       column_title_gp = gpar(fontsize = 16, fontface = "bold"),
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
  
  # Save PNG version
  png(file.path(figures_dir, "hierarchical_mp_clustering_heatmap.png"), 
      width = 16, height = 10, units = "in", res = 300)
  draw(combined_heatmap,
       column_title = "Meta-Program Hierarchical Clustering Analysis", 
       column_title_gp = gpar(fontsize = 16, fontface = "bold"),
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
  
  # Create a simpler version for presentations
  simple_heatmap <- Heatmap(
    t(patient_mp_matrix),
    name = "MP Proportion",
    cluster_rows = mp_hclust,
    cluster_columns = patient_hclust,
    col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    top_annotation = HeatmapAnnotation(
      MSI_Status = patient_annotations$Microsatellite_Status,
      Response = patient_annotations$Response,
      col = list(
        MSI_Status = msi_colors[names(msi_colors) %in% unique(patient_annotations$Microsatellite_Status)],
        Response = response_colors[names(response_colors) %in% unique(patient_annotations$Response)]
      )
    ),
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 12),
    column_title = "Patients",
    row_title = "Meta-Programs"
  )
  
  pdf(file.path(figures_dir, "hierarchical_mp_clustering_simple.pdf"), 
      width = 12, height = 8)
  draw(simple_heatmap,
       column_title = "Meta-Program Patient Clustering",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
  
  # ========================================================================
  # 9. SAVE CLUSTERING RESULTS AND STATISTICS
  # ========================================================================
  
  cat("  📋 Saving clustering results...\n")
  
  # Save Jaccard similarity matrix
  write.csv(jaccard_matrix, file.path(output_dir, "mp_jaccard_similarity_matrix.csv"))
  
  # Save patient clustering assignments
  patient_clusters <- cutree(patient_hclust, k = 4)  # Cut into 4 clusters
  mp_clusters <- cutree(mp_hclust, k = 3)  # Cut MPs into 3 clusters
  
  patient_cluster_df <- data.frame(
    Patient_ID = names(patient_clusters),
    Cluster = patient_clusters,
    stringsAsFactors = FALSE
  )
  write.csv(patient_cluster_df, file.path(output_dir, "patient_clusters.csv"), row.names = FALSE)
  
  mp_cluster_df <- data.frame(
    Meta_Program = names(mp_clusters),
    Cluster = mp_clusters,
    stringsAsFactors = FALSE
  )
  write.csv(mp_cluster_df, file.path(output_dir, "mp_clusters.csv"), row.names = FALSE)
  
  # Calculate clustering statistics
  clustering_stats <- list(
    n_patients = nrow(patient_mp_matrix),
    n_mps = ncol(patient_mp_matrix),
    avg_jaccard_similarity = mean(jaccard_matrix[upper.tri(jaccard_matrix)]),
    max_jaccard_similarity = max(jaccard_matrix[jaccard_matrix != 1]),
    min_jaccard_similarity = min(jaccard_matrix[jaccard_matrix != 1]),
    patient_cluster_sizes = table(patient_clusters),
    mp_cluster_sizes = table(mp_clusters)
  )
  
  # Save statistics
  capture.output({
    cat("Hierarchical Clustering Statistics\n")
    cat("=================================\n")
    cat("Number of patients:", clustering_stats$n_patients, "\n")
    cat("Number of meta-programs:", clustering_stats$n_mps, "\n")
    cat("Average Jaccard similarity:", round(clustering_stats$avg_jaccard_similarity, 3), "\n")
    cat("Max Jaccard similarity:", round(clustering_stats$max_jaccard_similarity, 3), "\n")
    cat("Min Jaccard similarity:", round(clustering_stats$min_jaccard_similarity, 3), "\n")
    cat("\nPatient cluster sizes:\n")
    print(clustering_stats$patient_cluster_sizes)
    cat("\nMP cluster sizes:\n")
    print(clustering_stats$mp_cluster_sizes)
  }, file = file.path(output_dir, "clustering_statistics.txt"))
  
  cat("  ✅ Hierarchical clustering heatmap completed!\n")
  cat("     📁 Files saved:\n")
  cat("        - hierarchical_mp_clustering_heatmap.pdf\n")
  cat("        - hierarchical_mp_clustering_heatmap.png\n") 
  cat("        - hierarchical_mp_clustering_simple.pdf\n")
  cat("        - mp_jaccard_similarity_matrix.csv\n")
  cat("        - patient_clusters.csv\n")
  cat("        - mp_clusters.csv\n")
  cat("        - clustering_statistics.txt\n")
  
  return(list(
    jaccard_matrix = jaccard_matrix,
    patient_clusters = patient_clusters,
    mp_clusters = mp_clusters,
    patient_hclust = patient_hclust,
    mp_hclust = mp_hclust,
    statistics = clustering_stats
  ))
}

# =============================================================================
# USAGE EXAMPLE
# =============================================================================

# Add this function to your main R analysis pipeline:
run_hierarchical_clustering_analysis <- function(data_dir, output_dir = NULL, figures_dir = NULL) {
  
  if (is.null(output_dir)) output_dir <- file.path(data_dir, "clustering_analysis")
  if (is.null(figures_dir)) figures_dir <- file.path(output_dir, "figures")
  
  
  # Load data
  cat("📥 Loading data for hierarchical clustering analysis...\n")
  patient_data <- fread(file.path(data_dir, "patient_level_summary.csv"))
  mp_signatures <- fread(file.path(data_dir, "mp_signature_genes.csv"))
  
  # Run hierarchical clustering analysis
  clustering_results <- create_hierarchical_mp_heatmap(patient_data, mp_signatures, output_dir, figures_dir)
  
  cat("🎉 Hierarchical clustering analysis complete!\n")
  return(clustering_results)
}

cat("📦 Hierarchical Clustering Heatmap Functions Loaded!\n")
cat("="*60, "\n")
cat("🔧 Available Functions:\n")
cat("  - create_hierarchical_mp_heatmap()\n")
cat("  - run_hierarchical_clustering_analysis()\n")
cat("\n💡 Quick Usage:\n")
cat("# Load your data\n")
cat("patient_data <- fread('patient_level_summary.csv')\n")
cat("mp_signatures <- fread('mp_signature_genes.csv')\n")
cat("\n# Create the heatmap\n")
cat("results <- create_hierarchical_mp_heatmap(patient_data, mp_signatures, output_dir, figures_dir)\n")
cat("\n# Or run complete analysis\n")
cat("results <- run_hierarchical_clustering_analysis(data_dir)\n")
cat("\n✨ Creates Figure 1D style hierarchical clustering heatmap!\n")