# R Functions for Loading Sparse Matrix Data from Python Export
# Optimized for single-cell data analysis
# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(rstatix)
library(broom)
library(effsize)
library(patchwork)
library(scales)
library(viridis)
library(ggalluvial)
library(gridExtra)
library(Matrix)
library(cluster)
library(dendextend)
library(ggrepel)

# Helper function to save plots with consistent settings
save_plot <- function(plot, figurePath, filename, width = 10, height = 8, dpi = 300) {
  ggsave(
    filename = file.path(figurePath, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "pdf"
  )

  print(paste("✅ Saved at:", file.path(figurePath, filename)))
}

# Helper function to add significance annotations
add_significance <- function(p_value) {
  case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# =============================================================================
# 1. DATA LOADING AND PREPROCESSING
# =============================================================================
load_sparse_expression_data <- function(data_dir, load_expression = TRUE) {
  cat("📥 Loading sparse matrix data from Python export...\n")
  cat("📁 Data directory:", data_dir, "\n\n")

  # Initialize results list
  data_list <- list()

  # 1. Load clinical metadata
  cat("👥 Loading clinical metadata...\n")
  if (file.exists(file.path(data_dir, "clinical_metadata.csv"))) {
    clinical_data <- fread(file.path(data_dir, "clinical_metadata.csv"))
    data_list$clinical_metadata <- clinical_data
    cat("   ✅ Clinical metadata:", nrow(clinical_data), "cells\n")
  }

  # 2. Load patient-level summary
  cat("📊 Loading patient-level summary...\n")
  if (file.exists(file.path(data_dir, "patient_level_summary.csv"))) {
    patient_data <- fread(file.path(data_dir, "patient_level_summary.csv"))
    data_list$patient_summary <- patient_data
    cat("   ✅ Patient summary:", nrow(patient_data), "patients\n")
  }

  # 3. Load MP signature genes
  cat("🧬 Loading MP signature genes...\n")
  if (file.exists(file.path(data_dir, "mp_signature_genes.csv"))) {
    mp_signatures <- fread(file.path(data_dir, "mp_signature_genes.csv"))
    data_list$mp_signatures <- mp_signatures
    cat("   ✅ MP signatures:", ncol(mp_signatures) - 1, "meta-programs\n")
  }

  # 4. Load sparse expression matrix (if requested)
  if (load_expression) {
    cat("🔬 Loading sparse expression matrix...\n")

    # Check if all required files exist
    mtx_file <- file.path(data_dir, "expression_matrix.mtx")
    genes_file <- file.path(data_dir, "gene_names.csv")
    barcodes_file <- file.path(data_dir, "cell_barcodes.csv")

    if (file.exists(mtx_file) && file.exists(genes_file) && file.exists(barcodes_file)) {
      # Load sparse matrix
      expression_matrix <- readMM(mtx_file)

      # Load gene names and cell barcodes
      gene_names <- fread(genes_file, header = TRUE)$gene_name
      cell_barcodes <- fread(barcodes_file, header = TRUE)$barcode

      # Set row and column names
      rownames(expression_matrix) <- cell_barcodes
      colnames(expression_matrix) <- gene_names

      # Convert to dgCMatrix for efficiency
      expression_matrix <- as(expression_matrix, "dgCMatrix")

      data_list$expression_matrix <- expression_matrix

      cat("   ✅ Expression matrix:", nrow(expression_matrix), "cells ×", ncol(expression_matrix), "genes\n")
      cat("   💾 Matrix sparsity:", round((1 - nnzero(expression_matrix) / length(expression_matrix)) * 100, 1), "%\n")

      # Load gene metadata if available
      if (file.exists(file.path(data_dir, "gene_metadata.csv"))) {
        gene_metadata <- fread(file.path(data_dir, "gene_metadata.csv"))
        data_list$gene_metadata <- gene_metadata
        cat("   ✅ Gene metadata loaded\n")
      }
    } else {
      cat("   ⚠️ Sparse matrix files not found, skipping expression data\n")
    }
  }

  # 5. Load MSI treated data
  cat("🎯 Loading MSI treated patient data...\n")
  if (file.exists(file.path(data_dir, "msi_treated_patients.csv"))) {
    msi_patients <- fread(file.path(data_dir, "msi_treated_patients.csv"))
    data_list$msi_patients <- msi_patients
    cat("   ✅ MSI treated patients:", nrow(msi_patients), "patients\n")
  }

  # 6. Load analysis summary
  if (file.exists(file.path(data_dir, "analysis_summary.csv"))) {
    analysis_summary <- fread(file.path(data_dir, "analysis_summary.csv"))
    data_list$analysis_summary <- analysis_summary
  }

  cat("\n📋 Data loading complete!\n")
  cat("Available datasets:\n")
  for (name in names(data_list)) {
    cat("  -", name, "\n")
  }

  # NEW: Load MP fraction data for stacked barplot
  cat("📊 Loading MP fraction data...\n")
  if (file.exists(file.path(data_dir, "sample_mp_fractions.csv"))) {
    mp_fractions <- fread(file.path(data_dir, "sample_mp_fractions.csv"))
    data_list$sample_mp_fractions <- mp_fractions
    cat("   ✅ MP fractions:", nrow(mp_fractions), "samples\n")
  }

  # NEW: Load MP marker genes
  cat("🧬 Loading MP marker genes...\n")
  if (file.exists(file.path(data_dir, "mp_marker_genes_ranked.csv"))) {
    mp_markers <- fread(file.path(data_dir, "mp_marker_genes_ranked.csv"))
    data_list$mp_marker_genes <- mp_markers
    cat("   ✅ MP markers:", nrow(mp_markers), "gene-MP pairs\n")
  }

  # NEW: Load expression subset for heatmap
  cat("🔥 Loading expression subset...\n")
  if (file.exists(file.path(data_dir, "expression_matrix_subset.csv"))) {
    expr_subset <- fread(file.path(data_dir, "expression_matrix_subset.csv"))
    # Convert first column to row names
    expr_matrix <- as.matrix(expr_subset[, -1])
    rownames(expr_matrix) <- expr_subset[[1]]
    data_list$expression_subset <- expr_matrix
    cat("   ✅ Expression subset:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "cells\n")
  }

  return(data_list)
}

# Create expression matrix subset for specific genes
subset_expression_matrix <- function(expression_matrix, genes) {
  # Find available genes
  available_genes <- intersect(genes, colnames(expression_matrix))
  missing_genes <- setdiff(genes, colnames(expression_matrix))

  if (length(missing_genes) > 0) {
    cat("⚠️ Missing genes:", length(missing_genes), "/", length(genes), "\n")
  }

  if (length(available_genes) == 0) {
    stop("No genes found in expression matrix")
  }

  # Subset matrix
  subset_matrix <- expression_matrix[, available_genes, drop = FALSE]

  cat("✅ Subsetted to", ncol(subset_matrix), "genes\n")
  return(subset_matrix)
}

# Calculate MP scores from expression matrix
calculate_mp_scores_from_expression <- function(expression_matrix, mp_signatures, method = "mean") {
  cat("🧮 Calculating MP scores from expression data...\n")

  # Get MP names (exclude gene column)
  mp_names <- setdiff(names(mp_signatures), c("V1", "gene", "Gene"))

  # Initialize results dataframe
  mp_scores <- data.frame(
    cell_barcode = rownames(expression_matrix),
    stringsAsFactors = FALSE
  )

  # Calculate scores for each MP
  for (mp_name in mp_names) {
    # Get signature genes for this MP
    signature_genes <- mp_signatures %>%
      filter(.data[[mp_name]] == 1) %>%
      pull(if ("gene" %in% names(.)) gene else V1)

    if (length(signature_genes) == 0) {
      cat("   ⚠️ No genes found for", mp_name, "\n")
      mp_scores[[paste0(mp_name, "_score")]] <- NA
      next
    }

    # Find available genes in expression matrix
    available_genes <- intersect(signature_genes, colnames(expression_matrix))

    if (length(available_genes) == 0) {
      cat("   ⚠️ No signature genes found in expression matrix for", mp_name, "\n")
      mp_scores[[paste0(mp_name, "_score")]] <- NA
      next
    }

    # Calculate scores
    if (length(available_genes) == 1) {
      scores <- expression_matrix[, available_genes]
    } else {
      if (method == "mean") {
        scores <- Matrix::rowMeans(expression_matrix[, available_genes])
      } else if (method == "median") {
        scores <- apply(expression_matrix[, available_genes], 1, median)
      } else if (method == "sum") {
        scores <- Matrix::rowSums(expression_matrix[, available_genes])
      } else {
        stop("Method must be 'mean', 'median', or 'sum'")
      }
    }

    mp_scores[[paste0(mp_name, "_score")]] <- scores
    cat("   ✅", mp_name, ":", length(available_genes), "genes used\n")
  }

  cat("✅ MP scores calculated for", length(mp_names), "meta-programs\n")
  return(mp_scores)
}

# Convert sparse matrix to dense for specific analyses
sparse_to_dense_safe <- function(sparse_matrix, max_genes = 2000) {
  n_genes <- ncol(sparse_matrix)
  n_cells <- nrow(sparse_matrix)

  # Estimate memory usage (rough)
  estimated_memory_gb <- (n_cells * n_genes * 8) / (1024^3) # 8 bytes per double

  cat("Matrix dimensions:", n_cells, "cells ×", n_genes, "genes\n")
  cat("Estimated memory for dense matrix:", round(estimated_memory_gb, 2), "GB\n")

  if (n_genes > max_genes) {
    cat("⚠️ Matrix too large for safe conversion to dense format\n")
    cat("Consider subsetting to <", max_genes, "genes first\n")
    return(NULL)
  }

  if (estimated_memory_gb > 4) {
    cat("⚠️ Large memory requirement detected\n")
    response <- readline("Continue with conversion? (y/n): ")
    if (tolower(response) != "y") {
      return(NULL)
    }
  }

  cat("Converting to dense matrix...\n")
  dense_matrix <- as.matrix(sparse_matrix)
  cat("✅ Conversion complete\n")

  return(dense_matrix)
}

# =============================================================================
# 2. Hierarchical Clustering Heatmap (Figure 1D Style)
# =============================================================================

create_marker_gene_correlation_heatmap <- function(data_dir, max_cells = 200,
                                                   top_genes_per_mp = 50,
                                                   output_dir = NULL, figures_dir = NULL) {
  cat("🔬 Creating Ordered MP Correlation Heatmap (MP1-MP10)...\n")
  cat("   Purpose: Validate MP assignments with ordered layout\n")

  # ========================================================================
  # 2. EXTRACT MARKER GENES FOR EACH MP
  # ========================================================================

  cat("  🧬 Extracting marker genes for each MP...\n")

  # Get MP columns (exclude gene name column)
  mp_cols <- names(mp_signatures)[!names(mp_signatures) %in% c("V1", "gene", "Gene")]

  if (length(mp_cols) == 0) {
    stop("❌ No valid MP columns found in mp_signatures")
  }

  # Extract marker genes for each MP
  mp_marker_genes <- list()

  for (mp_col in mp_cols) {
    # Get genes where this MP has value 1
    marker_genes <- mp_signatures %>%
      filter(.data[[mp_col]] == 1) %>%
      pull(if ("gene" %in% names(.)) gene else V1)

    # Clean MP name
    mp_name <- if (grepl("^[0-9]+$", mp_col)) paste0("MP", mp_col) else mp_col

    # Take top genes if too many
    if (length(marker_genes) > top_genes_per_mp) {
      marker_genes <- marker_genes[1:top_genes_per_mp]
    }

    mp_marker_genes[[mp_name]] <- marker_genes
    cat("     ", mp_name, ":", length(marker_genes), "marker genes\n")
  }

  # Get all unique marker genes
  all_marker_genes <- unique(unlist(mp_marker_genes))

  # Check which marker genes are available in expression matrix
  available_marker_genes <- intersect(all_marker_genes, colnames(expression_matrix))
  missing_genes <- setdiff(all_marker_genes, colnames(expression_matrix))

  cat("     Total unique marker genes:", length(all_marker_genes), "\n")
  cat("     Available in expression matrix:", length(available_marker_genes), "\n")
  if (length(missing_genes) > 0) {
    cat(
      "     Missing genes:", min(10, length(missing_genes)), "/", length(missing_genes),
      if (length(missing_genes) > 10) " (showing first 10)" else "", "\n"
    )
    cat("     ", paste(head(missing_genes, 10), collapse = ", "), "\n")
  }

  if (length(available_marker_genes) < 20) {
    stop("❌ Too few marker genes available in expression matrix (", length(available_marker_genes), ")")
  }

  # Update marker gene lists to only include available genes
  mp_marker_genes <- lapply(mp_marker_genes, function(genes) {
    intersect(genes, available_marker_genes)
  })

  # Remove MPs with too few marker genes
  mp_marker_genes <- mp_marker_genes[sapply(mp_marker_genes, length) >= 5]

  cat("     Final MPs with sufficient marker genes:", length(mp_marker_genes), "\n")
  for (mp_name in names(mp_marker_genes)) {
    cat("       ", mp_name, ":", length(mp_marker_genes[[mp_name]]), "genes\n")
  }

  # ========================================================================
  # 3. PREPARE CELL-LEVEL DATA
  # ========================================================================

  cat("  👥 Preparing cell-level data with MP assignments...\n")

  # Filter cells with valid MP assignments
  cells_with_mp <- clinical_data %>%
    filter(!is.na(MP_assignment) & MP_assignment != "Unresolved" & MP_assignment != "") %>%
    filter(!is.na(Patient_ID))

  # Match cell barcodes between clinical data and expression matrix
  common_cells <- intersect(rownames(cells_with_mp), rownames(expression_matrix))
  cells_with_mp <- cells_with_mp[match(common_cells, rownames(cells_with_mp)), ]

  cat("     Cells with MP assignments and expression data:", nrow(cells_with_mp), "\n")
  cat("     MP distribution:\n")
  print(table(cells_with_mp$MP_assignment))

  # Stratified downsampling by MP
  if (nrow(cells_with_mp) > max_cells) {
    cat("  📉 Downsampling to", max_cells, "cells (stratified by MP)...\n")

    # Calculate target cells per MP (proportional to current distribution)
    mp_counts <- table(cells_with_mp$MP_assignment)
    mp_props <- mp_counts / sum(mp_counts)
    target_per_mp <- pmax(50, round(max_cells * mp_props)) # At least 50 cells per MP

    # Adjust if total exceeds max_cells
    if (sum(target_per_mp) > max_cells) {
      target_per_mp <- round(target_per_mp * max_cells / sum(target_per_mp))
    }

    # Sample cells
    idx <- sample.int(nrow(cells_with_mp), size = max_cells, replace = FALSE)
    cells_with_mp <- cells_with_mp[idx, ]

    cat("     Final cell count:", nrow(cells_with_mp), "\n")
    cat("     Final MP distribution:\n")
    print(table(cells_with_mp$MP_assignment))
  }

  # ========================================================================
  # 4. ORDER CELLS BY MP ASSIGNMENT (MP1, MP2, ..., MP10)
  # ========================================================================

  cat("  📊 Ordering cells by MP assignment (MP1 to MP10)...\n")

  # Create ordered MP levels (MP1, MP2, ..., MP10)
  mp_levels <- paste0("MP", 1:10)
  available_mps <- intersect(mp_levels, unique(cells_with_mp$MP_assignment))

  # Order cells by MP assignment
  cells_with_mp$MP_assignment <- factor(cells_with_mp$MP_assignment, levels = available_mps)
  cells_with_mp <- cells_with_mp %>%
    arrange(MP_assignment)

  cat("     Cells ordered by MP assignment:", paste(available_mps, collapse = ", "), "\n")
  cat("     Final ordered MP distribution:\n")
  print(table(cells_with_mp$MP_assignment))

  # ========================================================================
  # 5. EXTRACT MARKER GENE EXPRESSION FOR ORDERED CELLS
  # ========================================================================

  cat("  🔬 Extracting marker gene expression for ordered cells...\n")

  # Get expression data for selected cells and marker genes
  selected_cell_barcodes <- rownames(cells_with_mp)
  marker_expression <- expression_matrix[selected_cell_barcodes, available_marker_genes]

  # Convert to dense matrix for correlation calculation
  cat("     Converting to dense matrix for correlation calculation...\n")
  marker_expression_dense <- as.matrix(marker_expression)

  cat("     Marker gene expression matrix:", dim(marker_expression_dense), "\n")
  cat("     Matrix sparsity:", round((1 - nnzero(marker_expression) / length(marker_expression)) * 100, 1), "%\n")

  # ========================================================================
  # 6. CALCULATE CELL-CELL CORRELATIONS
  # ========================================================================

  cat("  🧮 Calculating cell-cell correlations based on marker gene expression...\n")

  # Calculate correlation matrix between cells
  cat("     Computing correlations for", nrow(marker_expression_dense), "cells...\n")

  # Use Pearson correlation
  cell_correlation_matrix <- cor(t(marker_expression_dense))

  # Handle any NaN values
  cell_correlation_matrix[is.na(cell_correlation_matrix)] <- 0

  # Ensure row and column names match the ordered cells
  rownames(cell_correlation_matrix) <- selected_cell_barcodes
  colnames(cell_correlation_matrix) <- selected_cell_barcodes

  cat("  ✅ Cell-cell correlation matrix calculated\n")
  cat("     Matrix dimensions:", dim(cell_correlation_matrix), "\n")
  cat("     Average correlation:", round(mean(cell_correlation_matrix[upper.tri(cell_correlation_matrix)]), 3), "\n")

  # ========================================================================
  # 7. PREPARE ANNOTATIONS (NO CLUSTERING - KEEP MP ORDER)
  # ========================================================================

  cat("  🏷️ Preparing cell annotations...\n")

  # Create cell annotations matching the correlation matrix order
  cell_annotations <- cells_with_mp[, c("MP_assignment", "Treatment_Strategy", "Microsatellite_Status", "Response")]

  # Clean annotations
  cell_annotations <- cell_annotations %>%
    mutate(
      MP_assignment = factor(MP_assignment, levels = available_mps), # Keep MP order
      Treatment_Strategy = ifelse(is.na(Treatment_Strategy) | Treatment_Strategy == "",
        "untreated", Treatment_Strategy
      ),
      Microsatellite_Status = ifelse(is.na(Microsatellite_Status) | Microsatellite_Status == "",
        "Unknown", Microsatellite_Status
      ),
      Response = ifelse(is.na(Response) | Response == "", "Unknown", Response)
    )

  # ========================================================================
  # 8. CREATE COLOR SCHEMES
  # ========================================================================

  cat("  🎨 Setting up color schemes...\n")

  # MP colors - use distinct colors for each MP in order
  mp_colors <- setNames(
    rainbow(length(available_mps), start = 0, end = 0.85), # Avoid red at the end
    available_mps
  )

  # Clinical colors
  treatment_colors <- c(
    "Anti-PD1" = "#E74C3C",
    "Anti-PD1 plus Celecoxib" = "#3498DB",
    "Anti-PD1 plus CapeOx" = "#F39C12",
    "untreated" = "#95A5A6"
  )

  msi_colors <- c("MSI" = "#E74C3C", "MSS" = "#3498DB", "Unknown" = "#95A5A6")
  response_colors <- c("pCR" = "#27AE60", "non_pCR" = "#E67E22", "Unknown" = "#95A5A6")

  # ========================================================================
  # 9. CREATE ANNOTATIONS
  # ========================================================================

  # Row annotation (left) - MP assignments
  row_annotation <- rowAnnotation(
    MP_Assignment = cell_annotations$MP_assignment,
    col = list(
      MP_Assignment = mp_colors
    ),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      MP_Assignment = list(
        title = "Meta-Program",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        grid_height = unit(4, "mm")
      )
    ),
    annotation_width = unit(1.5, "cm")
  )

  # Column annotation (top) - clinical variables
  col_annotation <- HeatmapAnnotation(
    Treatment = cell_annotations$Treatment_Strategy,
    MSI_Status = cell_annotations$Microsatellite_Status,
    Response = cell_annotations$Response,
    col = list(
      Treatment = treatment_colors[names(treatment_colors) %in% unique(cell_annotations$Treatment_Strategy)],
      MSI_Status = msi_colors[names(msi_colors) %in% unique(cell_annotations$Microsatellite_Status)],
      Response = response_colors[names(response_colors) %in% unique(cell_annotations$Response)]
    ),
    annotation_name_gp = gpar(fontsize = 9),
    annotation_height = unit(1, "cm"),
    annotation_legend_param = list(
      Treatment = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)),
      MSI_Status = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)),
      Response = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8))
    )
  )

  # ========================================================================
  # 10. CREATE MAIN HEATMAP (NO CLUSTERING, PURPLE-YELLOW COLORS)
  # ========================================================================

  cat("  🔥 Creating ordered MP correlation heatmap (purple-yellow)...\n")

  # Create the main heatmap
  main_heatmap <- Heatmap(
    cell_correlation_matrix,
    name = "Cell-Cell\nCorrelation\n(Marker Genes)",

    # NO CLUSTERING - keep MP order
    cluster_rows = FALSE,
    cluster_columns = FALSE,

    # Purple to Yellow color scheme (0 to 1)
    col = colorRamp2(c(0, 0.5, 1), c("#440154", "#FDE725", "#FDE725")), # Purple to Yellow

    # Annotations
    left_annotation = row_annotation,
    top_annotation = col_annotation,

    # Appearance
    show_row_names = FALSE,
    show_column_names = FALSE,

    # Heatmap body
    rect_gp = gpar(col = "white", lwd = 0.1),

    # Titles
    column_title = "Cells (Ordered by MP: MP1 → MP10)",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title = "Cells (Ordered by MP: MP1 → MP10)",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),

    # Legend
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_direction = "vertical",
      legend_height = unit(6, "cm")
    )
  )

  # ========================================================================
  # 11. SAVE HEATMAPS
  # ========================================================================

  cat("  💾 Saving heatmaps...\n")

  # Main comprehensive heatmap
  pdf(file.path(figures_dir, "ordered_mp_correlation_heatmap.pdf"),
    width = 10, height = 8
  )
  draw(main_heatmap,
    column_title = "Ordered MP Correlation Analysis (Purple-Yellow)",
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  dev.off()

  return(NULL)
}

# =============================================================================
# 3. DESCRIPTIVE STATISTICS AND DATA EXPLORATION
# =============================================================================

generate_descriptive_stats <- function(clinical_data, patient_data, output_dir, figures_dir) {
  cat("📊 Generating descriptive statistics for your data structure...\n")

  # Create figures directory if it doesn't exist
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  # Overall summary statistics
  overall_summary <- list(
    n_cells = nrow(clinical_data),
    n_patients = length(unique(clinical_data$Patient_ID)),
    n_studies = length(unique(clinical_data$study)),
    n_samples = length(unique(clinical_data$Sample_ID))
  )

  # Patient-level characteristics
  patient_summary <- patient_data %>%
    summarise(
      n_patients = n(),
      median_cells_per_patient = median(n_cells, na.rm = TRUE),
      mean_cells_per_patient = mean(n_cells, na.rm = TRUE),
      min_cells_per_patient = min(n_cells, na.rm = TRUE),
      max_cells_per_patient = max(n_cells, na.rm = TRUE),
      n_with_response = sum(!is.na(Response) & Response != "", na.rm = TRUE),
      n_msi = sum(Microsatellite_Status == "MSI", na.rm = TRUE),
      n_mss = sum(Microsatellite_Status == "MSS", na.rm = TRUE),
      n_pcr = sum(Response == "pCR", na.rm = TRUE),
      n_non_pcr = sum(Response == "non_pCR", na.rm = TRUE)
    )

  # Treatment characteristics
  treatment_summary <- patient_data %>%
    filter(!is.na(Treatment_Strategy) & Treatment_Strategy != "") %>%
    count(Treatment_Strategy, sort = TRUE)

  # Create comprehensive descriptive plots

  # 1. Patient Distribution Overview
  p1 <- patient_data %>%
    ggplot(aes(x = n_cells)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(aes(xintercept = median(n_cells, na.rm = TRUE)),
      color = "red", linetype = "dashed", size = 1
    ) +
    labs(
      title = "Distribution of Cells per Patient",
      subtitle = paste("Median:", round(median(patient_data$n_cells, na.rm = TRUE)), "cells"),
      x = "Number of Cells", y = "Number of Patients"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # 2. Treatment Strategy Distribution (Pie Chart)
  treatment_clean <- patient_data %>%
    filter(!is.na(Treatment_Strategy) & Treatment_Strategy != "") %>%
    count(Treatment_Strategy) %>%
    mutate(
      percentage = n / sum(n) * 100,
      label = paste0(Treatment_Strategy, "\n(", round(percentage, 1), "%)")
    )

  p2 <- ggplot(treatment_clean, aes(x = "", y = n, fill = Treatment_Strategy)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_viridis_d() +
    labs(title = "Treatment Strategy Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(title = "Treatment Strategy"))

  # 3. MSI Status Distribution (Bar Plot)
  msi_clean <- patient_data %>%
    filter(!is.na(Microsatellite_Status) & Microsatellite_Status != "") %>%
    count(Microsatellite_Status) %>%
    mutate(percentage = n / sum(n) * 100)

  p3 <- ggplot(msi_clean, aes(x = Microsatellite_Status, y = n, fill = Microsatellite_Status)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
      vjust = -0.5, size = 4
    ) +
    scale_fill_manual(values = c("MSI" = "#E74C3C", "MSS" = "#3498DB")) +
    labs(
      title = "Microsatellite Status Distribution",
      x = "Microsatellite Status", y = "Number of Patients"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  # 4. Treatment Response Distribution
  response_clean <- patient_data %>%
    filter(!is.na(Response) & Response != "" & Response %in% c("pCR", "non_pCR")) %>%
    count(Response) %>%
    mutate(percentage = n / sum(n) * 100)

  p4 <- ggplot(response_clean, aes(x = Response, y = n, fill = Response)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
      vjust = -0.5, size = 4
    ) +
    scale_fill_manual(values = c("pCR" = "#27AE60", "non_pCR" = "#E67E22")) +
    labs(
      title = "Treatment Response Distribution",
      x = "Treatment Response", y = "Number of Patients"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  # 5. Age Distribution by Gender
  age_gender_clean <- patient_data %>%
    filter(!is.na(Age) & !is.na(Gender) & Gender != "")

  p5 <- ggplot(age_gender_clean, aes(x = Gender, y = Age, fill = Gender)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("F" = "#E91E63", "M" = "#2196F3")) +
    labs(
      title = "Age Distribution by Gender",
      x = "Gender", y = "Age (years)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red")

  # 6. Study Composition
  study_counts <- patient_data %>%
    count(study) %>%
    mutate(percentage = n / sum(n) * 100)

  p6 <- ggplot(study_counts, aes(x = reorder(study, n), y = n)) +
    geom_bar(stat = "identity", fill = "darkorange", alpha = 0.8) +
    geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")),
      hjust = -0.1, size = 3.5
    ) +
    coord_flip() +
    labs(
      title = "Patient Distribution by Study",
      x = "Study", y = "Number of Patients"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # Combine plots
  combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6)
  combined_plot <- combined_plot + plot_annotation(
    title = "MSI-H CRC Dataset: Descriptive Statistics Overview",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

  # Save plots
  ggsave(file.path(figures_dir, "descriptive_statistics_overview.pdf"),
    combined_plot,
    width = 16, height = 20, dpi = 300
  )
  ggsave(file.path(figures_dir, "descriptive_statistics_overview.png"),
    combined_plot,
    width = 16, height = 20, dpi = 300
  )

  # Save summary statistics
  write_csv(
    data.frame(
      Metric = names(overall_summary),
      Value = unlist(overall_summary)
    ),
    file.path(output_dir, "overall_summary_stats.csv")
  )

  write_csv(patient_summary, file.path(output_dir, "patient_level_summary_stats.csv"))
  write_csv(treatment_summary, file.path(output_dir, "treatment_summary_stats.csv"))

  cat("✅ Descriptive statistics generated and saved\n")
  return(list(
    overall = overall_summary,
    patient = patient_summary,
    treatment = treatment_summary
  ))
}

# =============================================================================
# 4. META-PROGRAM DISTRIBUTION ANALYSIS
# =============================================================================

analyze_mp_distributions <- function(clinical_data, patient_data, figures_dir, output_dir) {
  cat("🔬 Analyzing Meta-Program distributions...\n")

  # Get MP assignment data
  mp_assignment_data <- as.data.frame(clinical_data) %>%
    filter(!is.na(MP_assignment) & MP_assignment != "" & MP_assignment != "Unresolved") %>%
    group_by(Patient_ID) %>%
    slice(1) %>% # One per patient
    ungroup()

  results <- list()

  # 1. MP Distribution by Treatment Strategy
  cat("  📈 MP distribution by Treatment Strategy...\n")

  treatment_mp_data <- mp_assignment_data %>%
    filter(!is.na(Treatment_Strategy) & Treatment_Strategy != "")

  if (nrow(treatment_mp_data) > 0) {
    # Chi-square test
    chi_test_treatment <- chisq.test(table(treatment_mp_data$Treatment_Strategy, treatment_mp_data$MP_assignment))

    # Effect size (Cramér's V)
    cramers_v_treatment <- sqrt(chi_test_treatment$statistic /
      (nrow(treatment_mp_data) *
        (min(
          length(unique(treatment_mp_data$Treatment_Strategy)),
          length(unique(treatment_mp_data$MP_assignment))
        ) - 1)))

    results$treatment_response <- list(
      chi_square = chi_test_treatment,
      cramers_v = cramers_v_treatment,
      crosstab = table(treatment_mp_data$Treatment_Strategy, treatment_mp_data$MP_assignment)
    )

    # Create stacked bar plot
    treatment_mp_summary <- treatment_mp_data %>%
      count(Treatment_Strategy, MP_assignment) %>%
      group_by(Treatment_Strategy) %>%
      mutate(percentage = n / sum(n) * 100)

    p1 <- ggplot(treatment_mp_summary, aes(x = Treatment_Strategy, y = percentage, fill = MP_assignment)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_viridis_d(name = "Meta-Program") +
      labs(
        title = "Meta-Program Distribution by Treatment Strategy",
        subtitle = paste(
          "Chi-square p =", round(chi_test_treatment$p.value, 4),
          "; Cramér's V =", round(cramers_v_treatment, 3)
        ),
        x = "Treatment Strategy", y = "Percentage (%)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      geom_text(aes(label = ifelse(percentage > 5, paste0(round(percentage, 1), "%"), "")),
        position = position_stack(vjust = 0.5), size = 3
      )

    ggsave(file.path(figures_dir, "mp_distribution_by_treatment.pdf"), p1, width = 12, height = 8)
  }

  # 2. MP Distribution by Microsatellite Status
  cat("  🧬 MP distribution by Microsatellite Status...\n")

  msi_mp_data <- mp_assignment_data %>%
    filter(!is.na(Microsatellite_Status) & Microsatellite_Status != "")

  if (nrow(msi_mp_data) > 0) {
    chi_test_msi <- chisq.test(table(msi_mp_data$Microsatellite_Status, msi_mp_data$MP_assignment))
    cramers_v_msi <- sqrt(chi_test_msi$statistic /
      (nrow(msi_mp_data) *
        (min(
          length(unique(msi_mp_data$Microsatellite_Status)),
          length(unique(msi_mp_data$MP_assignment))
        ) - 1)))

    results$msi_status <- list(
      chi_square = chi_test_msi,
      cramers_v = cramers_v_msi,
      crosstab = table(msi_mp_data$Microsatellite_Status, msi_mp_data$MP_assignment)
    )

    # Create grouped bar plot
    msi_mp_summary <- msi_mp_data %>%
      count(Microsatellite_Status, MP_assignment) %>%
      group_by(Microsatellite_Status) %>%
      mutate(percentage = n / sum(n) * 100)

    p2 <- ggplot(msi_mp_summary, aes(x = MP_assignment, y = percentage, fill = Microsatellite_Status)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("MSI" = "#E74C3C", "MSS" = "#3498DB")) +
      labs(
        title = "Meta-Program Distribution by Microsatellite Status",
        subtitle = paste(
          "Chi-square p =", round(chi_test_msi$p.value, 4),
          "; Cramér's V =", round(cramers_v_msi, 3)
        ),
        x = "Meta-Program", y = "Percentage (%)",
        fill = "Microsatellite Status"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      geom_text(aes(label = round(percentage, 1)),
        position = position_dodge(width = 0.9), vjust = -0.5, size = 3
      )

    ggsave(file.path(figures_dir, "mp_distribution_by_msi_status.pdf"), p2, width = 10, height = 6)
  }

  # 3. Overall MP Assignment Distribution (Pie Chart)
  mp_overall <- mp_assignment_data %>%
    count(MP_assignment) %>%
    mutate(
      percentage = n / sum(n) * 100,
      label = paste0(MP_assignment, "\n(", round(percentage, 1), "%)")
    )

  p3 <- ggplot(mp_overall, aes(x = "", y = n, fill = MP_assignment)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_viridis_d() +
    labs(title = "Overall Meta-Program Assignment Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(title = "Meta-Program"))

  ggsave(file.path(figures_dir, "mp_overall_distribution.pdf"), p3, width = 10, height = 8)

  # Save statistical results
  if (!is.null(results$treatment_response)) {
    write_csv(
      as.data.frame(results$treatment_response$crosstab),
      file.path(output_dir, "mp_treatment_crosstab.csv")
    )
  }

  if (!is.null(results$msi_status)) {
    write_csv(
      as.data.frame(results$msi_status$crosstab),
      file.path(output_dir, "mp_msi_crosstab.csv")
    )
  }

  return(results)
}

# =============================================================================
# 5. TREATMENT RESPONSE ANALYSIS (PRIMARY ANALYSIS)
# =============================================================================

analyze_msi_pre_treatment_response <- function(patient_data, figures_dir, output_dir) {
  cat("🎯 Analyzing treatment response associations...\n")

  # Get MP proportion columns
  mp_prop_cols <- names(patient_data)[grepl("_proportion$", names(patient_data))]

  # Focus on patients with MSI pre-treatment data
  response_data <- patient_data %>%
    filter(Treatment_Strategy %in% c("Anti-PD1", "Anti-PD1 plus CapeOx", "Anti-PD1 plus Celecoxib")) %>%
    filter(Microsatellite_Status == "MSI") %>%
    filter(Response %in% c("pCR", "non_pCR")) %>%
    filter(Treatment_Stage == "Pre")

  cat("  📊 Patients with response data:", nrow(response_data), "\n")

  results <- list()

  if (nrow(response_data) > 10 && length(mp_prop_cols) > 0) {
    # Test each MP proportion against response
    mp_response_tests <- map_dfr(mp_prop_cols, function(mp_col) {
      test_data <- response_data %>%
        filter(!is.na(.data[[mp_col]]))

      # t-test
      t_test <- t.test(test_data[[mp_col]] ~ test_data$Response)

      # Effect size (Cohen's d)
      effect_size <- cohen.d(test_data[[mp_col]], test_data$Response)

      # Wilcoxon test
      wilcox_test <- wilcox.test(test_data[[mp_col]] ~ test_data$Response)

      # Group statistics
      group_stats <- test_data %>%
        group_by(Response) %>%
        summarise(
          n = n(),
          mean = mean(.data[[mp_col]], na.rm = TRUE),
          median = median(.data[[mp_col]], na.rm = TRUE),
          sd = sd(.data[[mp_col]], na.rm = TRUE),
          .groups = "drop"
        )

      tibble(
        mp_program = str_remove(mp_col, "_proportion"),
        t_test_p = t_test$p.value,
        t_test_statistic = t_test$statistic,
        wilcox_p = wilcox_test$p.value,
        cohens_d = effect_size$estimate,
        effect_magnitude = effect_size$magnitude,
        n_pcr = ifelse("pCR" %in% group_stats$Response,
          group_stats$n[group_stats$Response == "pCR"], 0
        ),
        n_non_pcr = ifelse("non_pCR" %in% group_stats$Response,
          group_stats$n[group_stats$Response == "non_pCR"], 0
        ),
        mean_pcr = ifelse("pCR" %in% group_stats$Response,
          group_stats$mean[group_stats$Response == "pCR"], NA
        ),
        mean_non_pcr = ifelse("non_pCR" %in% group_stats$Response,
          group_stats$mean[group_stats$Response == "non_pCR"], NA
        ),
        median_pcr = ifelse("pCR" %in% group_stats$Response,
          group_stats$median[group_stats$Response == "pCR"], NA
        ),
        median_non_pcr = ifelse("non_pCR" %in% group_stats$Response,
          group_stats$median[group_stats$Response == "non_pCR"], NA
        )
      )
    }, .id = "test_id")

    # Multiple testing correction
    if (nrow(mp_response_tests) > 0) {
      mp_response_tests <- mp_response_tests %>%
        mutate(
          t_test_p_adj = p.adjust(t_test_p, method = "fdr"),
          wilcox_p_adj = p.adjust(wilcox_p, method = "fdr"),
          significant_t = t_test_p_adj < 0.05,
          significant_wilcox = wilcox_p_adj < 0.05
        ) %>%
        arrange(t_test_p)

      # Save results
      write_csv(mp_response_tests, file.path(output_dir, "mp_treatment_response_analysis.csv"))
      results$mp_response_tests <- mp_response_tests

      # Create comprehensive box plots for all MPs
      plot_data <- response_data %>%
        select(Patient_ID, Response, all_of(mp_prop_cols)) %>%
        pivot_longer(
          cols = all_of(mp_prop_cols),
          names_to = "MP_program",
          values_to = "Proportion"
        ) %>%
        mutate(MP_program = str_remove(MP_program, "_proportion"))

      # Box plot with statistical annotations
      p_response <- ggplot(plot_data, aes(x = Response, y = Proportion, fill = Response)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
        facet_wrap(~MP_program, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c("pCR" = "#27AE60", "non_pCR" = "#E67E22")) +
        theme_minimal() +
        theme(
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        stat_compare_means(method = "t.test", label = "p.format", size = 3) +
        labs(
          title = "Meta-Program Proportions by Treatment Response",
          subtitle = "Individual patient data with statistical comparisons",
          y = "Meta-Program Proportion", x = "Treatment Response",
          fill = "Response"
        )

      ggsave(file.path(figures_dir, "mp_response_association_boxplots.pdf"),
        p_response,
        width = 16, height = 12
      )

      # Volcano plot for effect sizes
      volcano_data <- mp_response_tests %>%
        mutate(
          neg_log_p = -log10(t_test_p_adj),
          significant = ifelse(significant_t & abs(cohens_d) > 0.2, "Significant", "Not Significant"),
          label = ifelse(significant_t & abs(cohens_d) > 0.2, mp_program, "")
        )

      p_volcano <- ggplot(volcano_data, aes(x = cohens_d, y = neg_log_p, color = significant)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "blue") +
        scale_color_manual(values = c("Significant" = "#E74C3C", "Not Significant" = "#95A5A6")) +
        geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
        labs(
          title = "Meta-Program Association with Treatment Response",
          subtitle = "Effect Size vs Statistical Significance",
          x = "Cohen's d (Effect Size)", y = "-log10(Adjusted P-value)",
          color = "Significance"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      ggsave(file.path(figures_dir, "mp_response_volcano_plot.pdf"),
        p_volcano,
        width = 10, height = 8
      )

      # Significant MPs summary
      significant_mps <- mp_response_tests %>%
        filter(significant_t | significant_wilcox)

      if (nrow(significant_mps) > 0) {
        cat("  🌟 Found", nrow(significant_mps), "significant MPs associated with response\n")
        print(significant_mps %>% select(mp_program, t_test_p_adj, cohens_d, effect_magnitude))
      } else {
        cat("  ℹ️ No significant MPs found associated with treatment response\n")
      }
    }
  } else {
    cat("  ⚠️ Insufficient data for treatment response analysis\n")
  }

  return(results)
}

analyze_msi_non_PCR_treatment_response <- function(patient_data, figures_dir, output_dir) {
  cat("🎯 Analyzing treatment response associations...\n")

  # Get MP proportion columns
  mp_prop_cols <- names(patient_data)[grepl("_proportion$", names(patient_data))]

  # Focus on patients with MSI pre-treatment data
  response_data <- patient_data %>%
    filter(Microsatellite_Status == "MSI") %>%
    filter(Treatment_Stage == "Pre") %>%
    filter(Response %in% c("pCR", "non_pCR")) %>%
    filter(Treatment_Strategy %in% c("Anti-PD1", "Anti-PD1 plus CapeOx", "Anti-PD1 plus Celecoxib"))

  cat("  📊 Patients with response data:", nrow(response_data), "\n")

  results <- list()

  if (nrow(response_data) > 10 && length(mp_prop_cols) > 0) {
    # Test each MP proportion against response
    mp_response_tests <- map_dfr(mp_prop_cols, function(mp_col) {
      test_data <- response_data %>%
        filter(!is.na(.data[[mp_col]]))

      # t-test
      t_test <- t.test(test_data[[mp_col]] ~ test_data$Response)

      # Effect size (Cohen's d)
      effect_size <- cohen.d(test_data[[mp_col]], test_data$Response)

      # Wilcoxon test
      wilcox_test <- wilcox.test(test_data[[mp_col]] ~ test_data$Response)

      # Group statistics
      group_stats <- test_data %>%
        group_by(Response) %>%
        summarise(
          n = n(),
          mean = mean(.data[[mp_col]], na.rm = TRUE),
          median = median(.data[[mp_col]], na.rm = TRUE),
          sd = sd(.data[[mp_col]], na.rm = TRUE),
          .groups = "drop"
        )

      tibble(
        mp_program = str_remove(mp_col, "_proportion"),
        t_test_p = t_test$p.value,
        t_test_statistic = t_test$statistic,
        wilcox_p = wilcox_test$p.value,
        cohens_d = effect_size$estimate,
        effect_magnitude = effect_size$magnitude,
        n_pcr = ifelse("pCR" %in% group_stats$Response,
          group_stats$n[group_stats$Response == "pCR"], 0
        ),
        n_non_pcr = ifelse("non_pCR" %in% group_stats$Response,
          group_stats$n[group_stats$Response == "non_pCR"], 0
        ),
        mean_pcr = ifelse("pCR" %in% group_stats$Response,
          group_stats$mean[group_stats$Response == "pCR"], NA
        ),
        mean_non_pcr = ifelse("non_pCR" %in% group_stats$Response,
          group_stats$mean[group_stats$Response == "non_pCR"], NA
        ),
        median_pcr = ifelse("pCR" %in% group_stats$Response,
          group_stats$median[group_stats$Response == "pCR"], NA
        ),
        median_non_pcr = ifelse("non_pCR" %in% group_stats$Response,
          group_stats$median[group_stats$Response == "non_pCR"], NA
        )
      )
    }, .id = "test_id")

    # Multiple testing correction
    if (nrow(mp_response_tests) > 0) {
      mp_response_tests <- mp_response_tests %>%
        mutate(
          t_test_p_adj = p.adjust(t_test_p, method = "fdr"),
          wilcox_p_adj = p.adjust(wilcox_p, method = "fdr"),
          significant_t = t_test_p_adj < 0.05,
          significant_wilcox = wilcox_p_adj < 0.05
        ) %>%
        arrange(t_test_p)

      # Save results
      write_csv(mp_response_tests, file.path(output_dir, "mp_treatment_response_analysis.csv"))
      results$mp_response_tests <- mp_response_tests

      # Create comprehensive box plots for all MPs
      plot_data <- response_data %>%
        select(Patient_ID, Response, all_of(mp_prop_cols)) %>%
        pivot_longer(
          cols = all_of(mp_prop_cols),
          names_to = "MP_program",
          values_to = "Proportion"
        ) %>%
        mutate(MP_program = str_remove(MP_program, "_proportion"))

      # Box plot with statistical annotations
      p_response <- ggplot(plot_data, aes(x = Response, y = Proportion, fill = Response)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
        facet_wrap(~MP_program, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c("pCR" = "#27AE60", "non_pCR" = "#E67E22")) +
        theme_minimal() +
        theme(
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        stat_compare_means(method = "t.test", label = "p.format", size = 3) +
        labs(
          title = "Meta-Program Proportions by Treatment Response",
          subtitle = "Individual patient data with statistical comparisons",
          y = "Meta-Program Proportion", x = "Treatment Response",
          fill = "Response"
        )

      ggsave(file.path(figures_dir, "mp_response_association_boxplots.pdf"),
        p_response,
        width = 16, height = 12
      )

      # Volcano plot for effect sizes
      volcano_data <- mp_response_tests %>%
        mutate(
          neg_log_p = -log10(t_test_p_adj),
          significant = ifelse(significant_t & abs(cohens_d) > 0.2, "Significant", "Not Significant"),
          label = ifelse(significant_t & abs(cohens_d) > 0.2, mp_program, "")
        )

      p_volcano <- ggplot(volcano_data, aes(x = cohens_d, y = neg_log_p, color = significant)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "blue") +
        scale_color_manual(values = c("Significant" = "#E74C3C", "Not Significant" = "#95A5A6")) +
        geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
        labs(
          title = "Meta-Program Association with Treatment Response",
          subtitle = "Effect Size vs Statistical Significance",
          x = "Cohen's d (Effect Size)", y = "-log10(Adjusted P-value)",
          color = "Significance"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

      ggsave(file.path(figures_dir, "mp_response_volcano_plot.pdf"),
        p_volcano,
        width = 10, height = 8
      )

      # Significant MPs summary
      significant_mps <- mp_response_tests %>%
        filter(significant_t | significant_wilcox)

      if (nrow(significant_mps) > 0) {
        cat("  🌟 Found", nrow(significant_mps), "significant MPs associated with response\n")
        print(significant_mps %>% select(mp_program, t_test_p_adj, cohens_d, effect_magnitude))
      } else {
        cat("  ℹ️ No significant MPs found associated with treatment response\n")
      }
    }
  } else {
    cat("  ⚠️ Insufficient data for treatment response analysis\n")
  }

  return(results)
}

# =============================================================================
# 6. MSI vs MSS COMPARATIVE ANALYSIS
# =============================================================================

analyze_msi_vs_mss <- function(patient_data, figures_dir, output_dir) {
  cat("🔬 Analyzing MSI vs MSS differences...\n")

  # Get MP proportion columns
  mp_prop_cols <- names(patient_data)[grepl("_proportion$", names(patient_data))]

  # Filter patients with MSI status
  msi_comparison_data <- patient_data %>%
    filter(Microsatellite_Status %in% c("MSI", "MSS"))

  cat("  📊 Patients for MSI/MSS comparison:", nrow(msi_comparison_data), "\n")

  results <- list()

  if (nrow(msi_comparison_data) > 10 && length(mp_prop_cols) > 0) {
    # Test each MP proportion between MSI and MSS
    msi_mss_tests <- map_dfr(mp_prop_cols, function(mp_col) {
      if (mp_col %in% names(msi_comparison_data)) {
        test_data <- msi_comparison_data %>%
          filter(!is.na(.data[[mp_col]]))

        if (nrow(test_data) > 10) {
          # t-test
          t_test <- t.test(test_data[[mp_col]] ~ test_data$Microsatellite_Status)

          # Effect size
          effect_size <- cohen.d(test_data[[mp_col]], test_data$Microsatellite_Status)

          # Wilcoxon test
          wilcox_test <- wilcox.test(test_data[[mp_col]] ~ test_data$Microsatellite_Status)

          # Group statistics
          group_stats <- test_data %>%
            group_by(Microsatellite_Status) %>%
            summarise(
              n = n(),
              mean = mean(.data[[mp_col]], na.rm = TRUE),
              median = median(.data[[mp_col]], na.rm = TRUE),
              sd = sd(.data[[mp_col]], na.rm = TRUE),
              .groups = "drop"
            )

          tibble(
            mp_program = str_remove(mp_col, "_proportion"),
            t_test_p = t_test$p.value,
            wilcox_p = wilcox_test$p.value,
            cohens_d = effect_size$estimate,
            effect_magnitude = effect_size$magnitude,
            n_msi = ifelse("MSI" %in% group_stats$Microsatellite_Status,
              group_stats$n[group_stats$Microsatellite_Status == "MSI"], 0
            ),
            n_mss = ifelse("MSS" %in% group_stats$Microsatellite_Status,
              group_stats$n[group_stats$Microsatellite_Status == "MSS"], 0
            ),
            mean_msi = ifelse("MSI" %in% group_stats$Microsatellite_Status,
              group_stats$mean[group_stats$Microsatellite_Status == "MSI"], NA
            ),
            mean_mss = ifelse("MSS" %in% group_stats$Microsatellite_Status,
              group_stats$mean[group_stats$Microsatellite_Status == "MSS"], NA
            )
          )
        }
      }
    })

    # Multiple testing correction
    if (nrow(msi_mss_tests) > 0) {
      msi_mss_tests <- msi_mss_tests %>%
        mutate(
          t_test_p_adj = p.adjust(t_test_p, method = "fdr"),
          wilcox_p_adj = p.adjust(wilcox_p, method = "fdr"),
          significant = t_test_p_adj < 0.05
        ) %>%
        arrange(t_test_p)

      write_csv(msi_mss_tests, file.path(output_dir, "msi_vs_mss_analysis.csv"))
      results$msi_mss_tests <- msi_mss_tests

      # Visualization: Box plots
      plot_data <- msi_comparison_data %>%
        select(Patient_ID, Microsatellite_Status, all_of(mp_prop_cols)) %>%
        pivot_longer(
          cols = all_of(mp_prop_cols),
          names_to = "MP_program",
          values_to = "Proportion"
        ) %>%
        mutate(MP_program = str_remove(MP_program, "_proportion"))

      p_msi_mss <- ggplot(plot_data, aes(x = Microsatellite_Status, y = Proportion, fill = Microsatellite_Status)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
        facet_wrap(~MP_program, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c("MSI" = "#E74C3C", "MSS" = "#3498DB")) +
        theme_minimal() +
        theme(
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        stat_compare_means(method = "t.test", label = "p.format", size = 3) +
        labs(
          title = "Meta-Program Differences between MSI and MSS",
          subtitle = "Patient-level proportions with statistical comparisons",
          y = "Meta-Program Proportion", x = "Microsatellite Status",
          fill = "MSI Status"
        )

      ggsave(file.path(figures_dir, "msi_vs_mss_comparison.pdf"),
        p_msi_mss,
        width = 16, height = 12
      )

      # Heatmap of significant differences
      significant_mps <- msi_mss_tests %>%
        filter(significant) %>%
        select(mp_program, mean_msi, mean_mss, cohens_d, t_test_p_adj)

      if (nrow(significant_mps) > 0) {
        heatmap_data <- significant_mps %>%
          select(mp_program, mean_msi, mean_mss) %>%
          pivot_longer(
            cols = c(mean_msi, mean_mss),
            names_to = "Group", values_to = "Mean_Proportion"
          ) %>%
          mutate(Group = ifelse(Group == "mean_msi", "MSI", "MSS"))

        p_heatmap <- ggplot(heatmap_data, aes(x = Group, y = mp_program, fill = Mean_Proportion)) +
          geom_tile(color = "white", size = 0.5) +
          scale_fill_gradient2(
            low = "blue", mid = "white", high = "red",
            midpoint = median(heatmap_data$Mean_Proportion, na.rm = TRUE),
            name = "Mean\nProportion"
          ) +
          labs(
            title = "Significant Meta-Program Differences: MSI vs MSS",
            x = "Microsatellite Status", y = "Meta-Program"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold")
          )

        ggsave(file.path(figures_dir, "msi_vs_mss_heatmap.pdf"),
          p_heatmap,
          width = 8, height = 6
        )
      }

      cat("  📊 MSI vs MSS analysis completed:", nrow(msi_mss_tests), "MPs tested\n")
      if (nrow(significant_mps) > 0) {
        cat("  🌟 Found", nrow(significant_mps), "significant differences\n")
      }
    }
  }

  return(results)
}

# =============================================================================
# 7. COMPREHENSIVE RESULTS SUMMARY AND REPORT
# =============================================================================

generate_comprehensive_report <- function(mp_distribution_results, response_results,
                                          msi_mss_results, patient_data, output_dir, figures_dir) {
  cat("📋 Generating comprehensive analysis report...\n")

  # Get MP proportion columns for summary
  mp_prop_cols <- names(patient_data)[grepl("_proportion$", names(patient_data))]

  # Create comprehensive summary report
  report <- list(
    analysis_date = Sys.Date(),
    data_summary = list(
      n_patients = nrow(patient_data),
      n_msi_patients = sum(patient_data$Microsatellite_Status == "MSI", na.rm = TRUE),
      n_mss_patients = sum(patient_data$Microsatellite_Status == "MSS", na.rm = TRUE),
      n_response_data = sum(!is.na(patient_data$Response) &
        patient_data$Response %in% c("pCR", "non_pCR")),
      n_pcr = sum(patient_data$Response == "pCR", na.rm = TRUE),
      n_non_pcr = sum(patient_data$Response == "non_pCR", na.rm = TRUE),
      n_meta_programs = length(mp_prop_cols),
      treatment_strategies = unique(patient_data$Treatment_Strategy[
        !is.na(patient_data$Treatment_Strategy) & patient_data$Treatment_Strategy != ""
      ])
    ),
    statistical_tests_performed = list(
      mp_distribution_by_treatment = !is.null(mp_distribution_results$treatment_response),
      mp_distribution_by_msi = !is.null(mp_distribution_results$msi_status),
      treatment_response_analysis = !is.null(response_results$mp_response_tests),
      msi_vs_mss_analysis = !is.null(msi_mss_results$msi_mss_tests)
    ),
    significant_findings = list()
  )

  # Add significant findings from treatment response analysis
  if (!is.null(response_results$mp_response_tests)) {
    significant_response <- response_results$mp_response_tests %>%
      filter(significant_t | significant_wilcox)

    if (nrow(significant_response) > 0) {
      report$significant_findings$treatment_response <- list(
        n_significant_mps = nrow(significant_response),
        mp_programs = significant_response$mp_program,
        min_p_value = min(significant_response$t_test_p_adj, na.rm = TRUE),
        max_effect_size = max(abs(significant_response$cohens_d), na.rm = TRUE),
        effect_directions = significant_response %>%
          select(mp_program, cohens_d, effect_magnitude) %>%
          mutate(direction = ifelse(cohens_d > 0, "Higher in pCR", "Higher in non-pCR"))
      )
    }
  }

  # Add significant findings from MSI vs MSS analysis
  if (!is.null(msi_mss_results$msi_mss_tests)) {
    significant_msi <- msi_mss_results$msi_mss_tests %>%
      filter(significant)

    if (nrow(significant_msi) > 0) {
      report$significant_findings$msi_vs_mss <- list(
        n_significant_mps = nrow(significant_msi),
        mp_programs = significant_msi$mp_program,
        min_p_value = min(significant_msi$t_test_p_adj, na.rm = TRUE),
        max_effect_size = max(abs(significant_msi$cohens_d), na.rm = TRUE),
        effect_directions = significant_msi %>%
          select(mp_program, cohens_d, effect_magnitude) %>%
          mutate(direction = ifelse(cohens_d > 0, "Higher in MSI", "Higher in MSS"))
      )
    }
  }

  # Create summary visualization dashboard
  create_summary_dashboard <- function() {
    # 1. MP Score Distribution Heatmap
    mp_score_cols <- names(patient_data)[grepl("MP[0-9]+_score$", names(patient_data))]

    if (length(mp_score_cols) > 0) {
      heatmap_data <- patient_data %>%
        select(Patient_ID, all_of(mp_score_cols)) %>%
        pivot_longer(
          cols = all_of(mp_score_cols),
          names_to = "MP", values_to = "Score"
        ) %>%
        mutate(MP = str_extract(MP, "MP[0-9]+"))

      p1 <- ggplot(heatmap_data, aes(x = MP, y = Score)) +
        geom_violin(aes(fill = MP), alpha = 0.7) +
        geom_boxplot(width = 0.1, outlier.alpha = 0.5) +
        scale_fill_viridis_d() +
        labs(
          title = "Meta-Program Score Distributions",
          x = "Meta-Program", y = "Score"
        ) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    }

    # 2. Treatment Response Summary
    if (!is.null(response_results$mp_response_tests)) {
      effect_size_plot_data <- response_results$mp_response_tests %>%
        mutate(
          significance = ifelse(significant_t, "Significant", "Not Significant"),
          abs_effect = abs(cohens_d)
        ) %>%
        arrange(desc(abs_effect))

      p2 <- ggplot(effect_size_plot_data, aes(
        x = reorder(mp_program, abs_effect),
        y = abs_effect, fill = significance
      )) +
        geom_bar(stat = "identity", alpha = 0.8) +
        scale_fill_manual(values = c("Significant" = "#E74C3C", "Not Significant" = "#95A5A6")) +
        coord_flip() +
        labs(
          title = "Treatment Response: Effect Sizes by Meta-Program",
          x = "Meta-Program", y = "Absolute Cohen's d",
          fill = "Significance"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }

    # 3. MSI vs MSS Summary
    if (!is.null(msi_mss_results$msi_mss_tests)) {
      msi_effect_plot_data <- msi_mss_results$msi_mss_tests %>%
        mutate(
          significance = ifelse(significant, "Significant", "Not Significant"),
          abs_effect = abs(cohens_d)
        ) %>%
        arrange(desc(abs_effect))

      p3 <- ggplot(msi_effect_plot_data, aes(
        x = reorder(mp_program, abs_effect),
        y = abs_effect, fill = significance
      )) +
        geom_bar(stat = "identity", alpha = 0.8) +
        scale_fill_manual(values = c("Significant" = "#3498DB", "Not Significant" = "#95A5A6")) +
        coord_flip() +
        labs(
          title = "MSI vs MSS: Effect Sizes by Meta-Program",
          x = "Meta-Program", y = "Absolute Cohen's d",
          fill = "Significance"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }

    # Combine plots if they exist
    plots_list <- list()
    if (exists("p1")) plots_list <- append(plots_list, list(p1))
    if (exists("p2")) plots_list <- append(plots_list, list(p2))
    if (exists("p3")) plots_list <- append(plots_list, list(p3))

    if (length(plots_list) > 0) {
      if (length(plots_list) == 3) {
        combined_dashboard <- p1 / (p2 + p3)
      } else if (length(plots_list) == 2) {
        combined_dashboard <- plots_list[[1]] + plots_list[[2]]
      } else {
        combined_dashboard <- plots_list[[1]]
      }

      combined_dashboard <- combined_dashboard + plot_annotation(
        title = "MSI-H CRC Meta-Program Analysis: Summary Dashboard",
        theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
      )

      ggsave(file.path(figures_dir, "analysis_summary_dashboard.pdf"),
        combined_dashboard,
        width = 16, height = 12, dpi = 300
      )
      ggsave(file.path(figures_dir, "analysis_summary_dashboard.png"),
        combined_dashboard,
        width = 16, height = 12, dpi = 300
      )
    }
  }

  # Create the dashboard
  create_summary_dashboard()

  # Create summary table
  summary_table <- tibble(
    Analysis = c(
      "MP Distribution by Treatment", "MP Distribution by MSI Status",
      "Treatment Response Association", "MSI vs MSS Comparison"
    ),
    Performed = c(
      !is.null(mp_distribution_results$treatment_response),
      !is.null(mp_distribution_results$msi_status),
      !is.null(response_results$mp_response_tests),
      !is.null(msi_mss_results$msi_mss_tests)
    ),
    Significant_Findings = c(
      if (!is.null(mp_distribution_results$treatment_response)) {
        mp_distribution_results$treatment_response$chi_square$p.value < 0.05
      } else {
        FALSE
      },
      if (!is.null(mp_distribution_results$msi_status)) {
        mp_distribution_results$msi_status$chi_square$p.value < 0.05
      } else {
        FALSE
      },
      if (!is.null(response_results$mp_response_tests)) {
        sum(response_results$mp_response_tests$significant_t |
          response_results$mp_response_tests$significant_wilcox, na.rm = TRUE) > 0
      } else {
        FALSE
      },
      if (!is.null(msi_mss_results$msi_mss_tests)) {
        sum(msi_mss_results$msi_mss_tests$significant, na.rm = TRUE) > 0
      } else {
        FALSE
      }
    ),
    N_Significant_MPs = c(
      if (!is.null(mp_distribution_results$treatment_response)) {
        ifelse(mp_distribution_results$treatment_response$chi_square$p.value < 0.05,
          "Overall significant", "Not significant"
        )
      } else {
        "Not tested"
      },
      if (!is.null(mp_distribution_results$msi_status)) {
        ifelse(mp_distribution_results$msi_status$chi_square$p.value < 0.05,
          "Overall significant", "Not significant"
        )
      } else {
        "Not tested"
      },
      if (!is.null(response_results$mp_response_tests)) {
        sum(response_results$mp_response_tests$significant_t |
          response_results$mp_response_tests$significant_wilcox, na.rm = TRUE)
      } else {
        0
      },
      if (!is.null(msi_mss_results$msi_mss_tests)) {
        sum(msi_mss_results$msi_mss_tests$significant, na.rm = TRUE)
      } else {
        0
      }
    )
  )

  write_csv(summary_table, file.path(output_dir, "analysis_summary_table.csv"))

  # Save detailed report as text
  report_text <- capture.output({
    cat("MSI-H CRC Meta-Program Analysis Report\n")
    cat("=====================================\n\n")
    cat("Analysis Date:", as.character(report$analysis_date), "\n\n")

    cat("DATA SUMMARY:\n")
    cat("- Total Patients:", report$data_summary$n_patients, "\n")
    cat("- MSI Patients:", report$data_summary$n_msi_patients, "\n")
    cat("- MSS Patients:", report$data_summary$n_mss_patients, "\n")
    cat("- Patients with Response Data:", report$data_summary$n_response_data, "\n")
    cat("- pCR Patients:", report$data_summary$n_pcr, "\n")
    cat("- non-pCR Patients:", report$data_summary$n_non_pcr, "\n")
    cat("- Meta-Programs Analyzed:", report$data_summary$n_meta_programs, "\n\n")

    cat("TREATMENT STRATEGIES:\n")
    for (strategy in report$data_summary$treatment_strategies) {
      cat("-", strategy, "\n")
    }
    cat("\n")

    if (!is.null(report$significant_findings$treatment_response)) {
      cat("SIGNIFICANT TREATMENT RESPONSE FINDINGS:\n")
      cat(
        "- Number of significant MPs:",
        report$significant_findings$treatment_response$n_significant_mps, "\n"
      )
      cat(
        "- Significant MPs:",
        paste(report$significant_findings$treatment_response$mp_programs, collapse = ", "), "\n"
      )
      cat(
        "- Minimum p-value:",
        round(report$significant_findings$treatment_response$min_p_value, 6), "\n"
      )
      cat(
        "- Maximum effect size:",
        round(report$significant_findings$treatment_response$max_effect_size, 3), "\n\n"
      )
    }

    if (!is.null(report$significant_findings$msi_vs_mss)) {
      cat("SIGNIFICANT MSI vs MSS FINDINGS:\n")
      cat(
        "- Number of significant MPs:",
        report$significant_findings$msi_vs_mss$n_significant_mps, "\n"
      )
      cat(
        "- Significant MPs:",
        paste(report$significant_findings$msi_vs_mss$mp_programs, collapse = ", "), "\n"
      )
      cat(
        "- Minimum p-value:",
        round(report$significant_findings$msi_vs_mss$min_p_value, 6), "\n"
      )
      cat(
        "- Maximum effect size:",
        round(report$significant_findings$msi_vs_mss$max_effect_size, 3), "\n\n"
      )
    }
  })

  writeLines(report_text, file.path(output_dir, "comprehensive_analysis_report.txt"))

  cat("✅ Comprehensive report generated and saved\n")
  cat("📁 Key outputs:\n")
  cat("   - analysis_summary_dashboard.pdf/png\n")
  cat("   - analysis_summary_table.csv\n")
  cat("   - comprehensive_analysis_report.txt\n")

  return(report)
}

# =============================================================================
# ADD NEW R VISUALIZATION FUNCTIONS
# =============================================================================

# Add after your existing R functions:
create_mp_composition_barplot_existing_style <- function(data_list, save_path, mp_colors = NULL) {
  # Use your existing colors object and ggplot style
  mp_data <- data_list$sample_mp_fractions

  # Define appealing color palette for MPs
  if (is.null(mp_colors)) {
    # Get number of unique MPs
    n_mps <- length(unique(mp_data$MP))

    # Create a sophisticated color palette
    if (n_mps <= 8) {
      mp_colors <- c(
        "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
        "#6A3D9A", "#B15928", "#A6CEE3", "#FDBF6F"
      )[1:n_mps]
    } else {
      # For more than 8 MPs, use a extended palette
      mp_colors <- c(
        "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
        "#6A3D9A", "#B15928", "#A6CEE3", "#FDBF6F",
        "#FB9A99", "#CAB2D6", "#FFFF99", "#B2DF8A"
      )
      if (n_mps > 12) {
        # Generate additional colors using rainbow
        extra_colors <- rainbow(n_mps - 12, start = 0.1, end = 0.9)
        mp_colors <- c(mp_colors, extra_colors)
      }
    }
    names(mp_colors) <- paste0("MP", sort(as.numeric(gsub("MP", "", unique(mp_data$MP)))))
  }

  # Calculate total fraction per sample for validation
  total_fractions <- mp_data %>%
    group_by(Sample_ID) %>%
    summarise(total_fraction = sum(MP_fraction), .groups = "drop")

  # Create "Unassigned" category for samples that don't sum to 1
  mp_data_complete <- mp_data %>%
    left_join(total_fractions, by = "Sample_ID") %>%
    group_by(Sample_ID, Response, Treatment_Stage) %>%
    do({
      current_data <- .
      if (current_data$total_fraction[1] < 0.999) { # Account for floating point precision
        unassigned_fraction <- 1 - current_data$total_fraction[1]
        unassigned_row <- current_data[1, ] %>%
          mutate(
            MP = "Unassigned",
            MP_fraction = unassigned_fraction,
            mp_cell_count = round(unassigned_fraction * total_cells[1])
          )
        rbind(current_data, unassigned_row)
      } else {
        current_data
      }
    }) %>%
    ungroup() %>%
    select(-total_fraction)

  # Add grey color for unassigned
  if ("Unassigned" %in% mp_data_complete$MP) {
    mp_colors <- c(mp_colors, "Unassigned" = "grey75")
  }

  # Order samples by treatment response, then by dominant MP
  sample_order <- mp_data_complete %>%
    group_by(Sample_ID, Response, Treatment_Stage) %>%
    slice_max(order_by = MP_fraction, n = 1) %>%
    arrange(Response, Treatment_Stage, desc(MP_fraction)) %>%
    pull(Sample_ID) %>%
    unique()

  # Reorder factor levels
  mp_data_complete$Sample_ID <- factor(mp_data_complete$Sample_ID, levels = sample_order)
  mp_data_complete$MP <- as.factor(mp_data_complete$MP)

  # Create treatment response labels for x-axis grouping
  sample_labels <- mp_data_complete %>%
    select(Sample_ID, Response, Treatment_Stage) %>%
    distinct() %>%
    mutate(
      label_color = case_when(
        Response == "pCR" ~ "#2E8B57", # Sea green for responders
        Response == "non_pCR" ~ "#CD5C5C", # Indian red for non-responders
        TRUE ~ "#666666"
      ),
      stage_symbol = case_when(
        Treatment_Stage == "Pre" ~ "●",
        Treatment_Stage == "Post" ~ "▲",
        TRUE ~ "■"
      )
    )

  # Main plot
  p <- ggplot(mp_data_complete, aes(x = Sample_ID, y = MP_fraction, fill = MP)) +
    # Background panel for better contrast
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "white", color = NA, inherit.aes = FALSE
    ) +

    # Main stacked bars
    geom_col(width = 0.85, color = "white", size = 0.2, alpha = 0.9) +

    # Add horizontal reference lines
    geom_hline(
      yintercept = c(0.25, 0.5, 0.75, 1.0),
      color = "grey90", linetype = "dotted", size = 0.3
    ) +

    # Color scale
    scale_fill_manual(values = mp_colors, name = "Meta-Program") +

    # Y-axis formatting
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.02),
      expand = c(0, 0)
    ) +

    # Sophisticated theme
    theme_minimal(base_size = 11) +
    theme(
      # Panel and background
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey95", size = 0.3),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),

      # Axes
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8,
        color = "grey30"
      ),
      axis.text.y = element_text(size = 9, color = "grey30"),
      axis.title = element_text(size = 11, color = "grey20", face = "bold"),
      axis.line.x = element_line(color = "grey60", size = 0.3),
      axis.ticks.x = element_line(color = "grey60", size = 0.3),
      axis.ticks.y = element_blank(),

      # Title and labels
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 20),
        color = "grey20"
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        color = "grey40",
        margin = margin(b = 15)
      ),

      # Legend
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold", color = "grey20"),
      legend.text = element_text(size = 9, color = "grey30"),
      legend.key.size = unit(0.8, "cm"),
      legend.key = element_rect(color = "white", size = 0.3),
      legend.margin = margin(t = 20),
      legend.box.spacing = unit(0.5, "cm"),

      # Margins
      plot.margin = margin(20, 20, 20, 20)
    ) +

    # Labels
    labs(
      title = "Meta-Program Composition Across Patient Samples",
      subtitle = "Fraction of malignant cells assigned to each meta-program",
      x = "Sample ID",
      y = "Meta-Program Fraction",
      caption = paste0(
        "Total samples: ", length(unique(mp_data_complete$Sample_ID)),
        " | Meta-programs: ", length(unique(mp_data_complete$MP[mp_data_complete$MP != "Unassigned"]))
      )
    ) +

    # Arrange legend in multiple rows if many MPs
    guides(
      fill = guide_legend(
        nrow = ceiling(length(mp_colors) / 6),
        byrow = TRUE,
        override.aes = list(color = "white", size = 0.3)
      )
    )

  # Add response status indicators at bottom
  # Extract sample info for annotations
  sample_annotations <- mp_data_complete %>%
    select(Sample_ID, Response, Treatment_Stage) %>%
    distinct() %>%
    mutate(
      x_pos = as.numeric(factor(Sample_ID, levels = sample_order)),
      response_color = case_when(
        Response == "pCR" ~ "#2E8B57",
        Response == "non_pCR" ~ "#CD5C5C",
        TRUE ~ "#999999"
      )
    )

  # Add response indicators
  p <- p +
    # Response status dots
    geom_point(
      data = sample_annotations,
      aes(x = x_pos, y = -0.04, color = I(response_color)),
      size = 2,
      shape = 16,
      inherit.aes = FALSE
    ) +
    # Extend y-axis slightly to accommodate indicators
    coord_cartesian(ylim = c(-0.06, 1.02), clip = "off")

  # Print summary statistics
  cat("\n📊 Meta-Program Composition Summary:\n")
  cat("Total samples:", length(unique(mp_data_complete$Sample_ID)), "\n")
  cat("Meta-programs identified:", length(unique(mp_data_complete$MP[mp_data_complete$MP != "Unassigned"])), "\n")

  # Save using your existing save_plot function
  save_plot(p, save_path, "mp_composition_stacked_barplot", width = 25, height = 8)
}

create_mp_composition_barplot_modern <- function(mp_data, save_path) {
  mp_data <- data_list$sample_mp_fractions

  # Prepare data with unassigned fractions
  mp_data_complete <- mp_data %>%
    group_by(Sample_ID) %>%
    mutate(
      total_assigned = sum(MP_fraction),
      unassigned = pmax(0, 1 - total_assigned)
    ) %>%
    ungroup()

  # Add unassigned rows
  if (any(mp_data_complete$unassigned > 0.001)) {
    unassigned_data <- mp_data_complete %>%
      filter(unassigned > 0.001) %>%
      select(Sample_ID, unassigned, Response, Treatment_Stage) %>%
      distinct() %>%
      mutate(
        MP = "Unassigned",
        MP_fraction = unassigned
      ) %>%
      select(Sample_ID, MP, MP_fraction, Response, Treatment_Stage)

    mp_data_plot <- bind_rows(
      mp_data %>% select(Sample_ID, MP, MP_fraction, Response, Treatment_Stage),
      unassigned_data
    )
  } else {
    mp_data_plot <- mp_data %>% select(Sample_ID, MP, MP_fraction, Response, Treatment_Stage)
  }

  # Modern color palette
  n_mps <- length(unique(mp_data_plot$MP[mp_data_plot$MP != "Unassigned"]))
  mp_colors <- plasma(n = n_mps, end = 0.9)
  names(mp_colors) <- paste0("MP", sort(as.numeric(gsub("MP", "", unique(mp_data_plot$MP[mp_data_plot$MP != "Unassigned"])))))

  if ("Unassigned" %in% mp_data_plot$MP) {
    mp_colors <- c(mp_colors, "Unassigned" = "grey75")
  }

  # Order samples
  sample_order <- mp_data_plot %>%
    filter(MP != "Unassigned") %>%
    group_by(Sample_ID, Response) %>%
    slice_max(order_by = MP_fraction, n = 1) %>%
    arrange(Response, desc(MP_fraction)) %>%
    pull(Sample_ID) %>%
    unique()

  mp_data_plot$Sample_ID <- factor(mp_data_plot$Sample_ID, levels = sample_order)

  # Modern minimal plot
  p <- ggplot(mp_data_plot, aes(x = Sample_ID, y = MP_fraction, fill = MP)) +
    geom_col(width = 0.9, alpha = 0.8) +
    scale_fill_manual(values = mp_colors, name = "Meta-Program") +
    scale_y_continuous(
      labels = percent_format(),
      expand = c(0, 0),
      limits = c(0, 1)
    ) +
    theme_void() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(
        size = 11,
        angle = 90,
        vjust = 0.5,
        margin = margin(r = 15)
      ),
      legend.position = "bottom",
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 20)
      ),
      panel.grid.major.y = element_line(color = "white", size = 0.5),
      plot.background = element_rect(fill = "grey98", color = NA),
      panel.background = element_rect(fill = "grey98", color = NA)
    ) +
    labs(
      title = "Meta-Program Landscape",
      y = "Fraction"
    )

  # Save modern version
  save_plot(p, save_path, "mp_composition_modern.pdf", width = 25, height = 8)

  return(p)
}

create_mp_features_heatmap <- function(data_list, save_path, top_genes_per_mp = 15) {
  cat("🎨 Creating MP Features Heatmap (Using Proper Data Structure)...\n")

  # Load required data
  clinical_data <- data_list$clinical_metadata
  patient_data <- data_list$patient_summary
  mp_marker_genes <- data_list$mp_marker_genes
  sample_mp_fractions <- data_list$sample_mp_fractions
  expression_subset <- data_list$expression_subset

  # Process
  mp_marker_genes$Meta_Program <- paste0("MP", mp_marker_genes$Meta_Program)

  cat(sprintf("📊 Data dimensions:\n"))
  cat(sprintf("   Clinical data: %d cells × %d variables\n", nrow(clinical_data), ncol(clinical_data)))
  cat(sprintf("   Sample MP fractions: %d samples × %d variables\n", nrow(sample_mp_fractions), ncol(sample_mp_fractions)))
  if (!is.null(mp_marker_genes)) {
    cat(sprintf("   MP marker genes: %d entries\n", nrow(mp_marker_genes)))
  }
  if (!is.null(expression_subset)) {
    cat(sprintf("   Expression subset: %d genes × %d cells\n", nrow(expression_subset), ncol(expression_subset)))
  }

  # =========================================================================
  # PART 1: Extract MP Information
  # =========================================================================

  # Get MP names from sample fractions (most reliable source)
  mp_cols <- grep("^MP[0-9]+$", unique(sample_mp_fractions$MP), value = TRUE)
  mp_names <- sort(mp_cols)

  cat(sprintf("   Found %d Meta-Programs: %s\n", length(mp_names), paste(mp_names, collapse = ", ")))

  # =========================================================================
  # PART 2: Prepare Gene Expression Matrix for Top Marker Genes
  # =========================================================================

  expression_heatmap_data <- NULL
  gene_annotation <- NULL

  if (!is.null(mp_marker_genes) && !is.null(expression_subset)) {
    cat("   Processing MP marker genes and expression data...\n")

    # Get top genes per MP from marker genes data
    top_genes_list <- list()

    if ("Meta_Program" %in% colnames(mp_marker_genes) && "gene" %in% colnames(mp_marker_genes)) {
      for (mp in mp_names) {
        mp_genes <- mp_marker_genes$gene[mp_marker_genes$Meta_Program == mp]

        # If we have ranking information, use it
        if ("rank" %in% colnames(mp_marker_genes)) {
          mp_gene_data <- mp_marker_genes[mp_marker_genes$Meta_Program == mp, ]
          mp_gene_data <- mp_gene_data[order(mp_gene_data$rank), ]
          top_genes_list[[mp]] <- head(mp_gene_data$gene, top_genes_per_mp)
        } else {
          top_genes_list[[mp]] <- head(mp_genes, top_genes_per_mp)
        }
      }
    }

    # Get all unique top genes
    all_top_genes <- unique(unlist(top_genes_list))

    # Filter to genes available in expression data
    available_genes <- intersect(all_top_genes, rownames(expression_subset))

    if (length(available_genes) > 5) {
      cat(sprintf(
        "   Using %d/%d marker genes with expression data\n",
        length(available_genes), length(all_top_genes)
      ))

      # Calculate mean expression per MP using sample fractions as weights
      expression_heatmap_data <- calculate_mp_weighted_expression(
        expression_subset, sample_mp_fractions, clinical_data,
        mp_names, available_genes
      )

      # Create gene annotation (which MP each gene primarily belongs to)
      gene_annotation <- create_gene_annotation_from_markers(
        available_genes, top_genes_list, mp_names
      )
    } else {
      cat("   ⚠️ Too few overlapping genes, skipping expression heatmap\n")
    }
  }

  # =========================================================================
  # PART 3: Calculate Clinical Associations from Sample Data
  # =========================================================================

  clinical_associations <- calculate_clinical_associations_from_samples(
    sample_mp_fractions, mp_names
  )

  # =========================================================================
  # PART 4: Calculate Sample-level Statistics
  # =========================================================================

  sample_stats <- calculate_sample_level_stats(sample_mp_fractions, mp_names)

  # =========================================================================
  # PART 5: Create ComplexHeatmap
  # =========================================================================

  # Color schemes
  clinical_colors <- circlize::colorRamp2(c(-0.5, 0, 0.5), c("#D73027", "white", "#1A9850"))
  expression_colors <- circlize::colorRamp2(c(-2, 0, 2), c("#053061", "white", "#67001F"))
  fraction_colors <- circlize::colorRamp2(c(0, 0.5, 1), c("white", "#FDCC8A", "#E31A1C"))

  # -------------------------------------------------------------------------
  # Panel 1: Clinical Associations (if available)
  # -------------------------------------------------------------------------

  clinical_heatmap <- NULL
  if (!is.null(clinical_associations) && nrow(clinical_associations) > 0) {
    clinical_heatmap <- Heatmap(
      as.matrix(clinical_associations),
      name = "Clinical\nAssociation",
      col = clinical_colors,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 10, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        value <- as.matrix(clinical_associations)[i, j]
        if (!is.na(value) && abs(value) > 0.02) {
          grid.text(sprintf("%.2f", value), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
        }
      },
      height = unit(nrow(clinical_associations) * 0.7, "cm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        grid_height = unit(0.4, "cm")
      )
    )
  }

  # -------------------------------------------------------------------------
  # Panel 2: Sample Statistics Annotations
  # -------------------------------------------------------------------------

  # Number of samples per MP
  n_samples_anno <- anno_barplot(
    sample_stats$n_samples,
    bar_width = 0.8,
    gp = gpar(fill = "#4292C6", col = "#08519C"),
    height = unit(1.2, "cm"),
    axis_param = list(gp = gpar(fontsize = 8))
  )

  # Mean MP fraction
  mean_fraction_anno <- anno_barplot(
    sample_stats$mean_fraction,
    bar_width = 0.8,
    gp = gpar(fill = "#74C476", col = "#238B45"),
    height = unit(1.0, "cm"),
    axis_param = list(gp = gpar(fontsize = 8))
  )

  # -------------------------------------------------------------------------
  # Panel 3: Main Heatmap
  # -------------------------------------------------------------------------

  if (!is.null(expression_heatmap_data)) {
    # Use expression data
    main_matrix <- expression_heatmap_data
    main_colors <- expression_colors
    main_name <- "Z-score\nExpression"

    # Z-score normalize across genes
    main_matrix_scaled <- t(scale(t(main_matrix)))
    main_matrix_scaled[is.infinite(main_matrix_scaled)] <- 0
    main_matrix_scaled[is.na(main_matrix_scaled)] <- 0

    # Create gene annotation colors
    gene_colors <- create_mp_color_scheme(mp_names)

    main_heatmap <- Heatmap(
      main_matrix_scaled,
      name = main_name,
      col = main_colors,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "right",
      column_names_side = "bottom",
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 11, fontface = "bold"),
      left_annotation = rowAnnotation(
        "Primary MP" = gene_annotation,
        col = list("Primary MP" = gene_colors),
        annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
        simple_anno_size = unit(0.4, "cm")
      ),
      top_annotation = HeatmapAnnotation(
        "# Samples" = n_samples_anno,
        "Mean Fraction" = mean_fraction_anno,
        annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
      ),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        grid_height = unit(0.4, "cm")
      )
    )
  } else {
    # Fallback: Use sample MP fractions as heatmap
    cat("   Creating sample fraction heatmap as main panel...\n")

    # Get sample fraction matrix
    sample_matrix <- as.matrix(sample_mp_fractions[, mp_names])
    rownames(sample_matrix) <- sample_mp_fractions$Sample_ID

    # Take subset for visualization if too many samples
    if (nrow(sample_matrix) > 50) {
      # Select samples with highest diversity (entropy)
      sample_entropy <- apply(sample_matrix, 1, function(x) {
        x <- x[x > 0]
        if (length(x) > 1) -sum(x * log(x)) else 0
      })
      top_samples <- head(order(sample_entropy, decreasing = TRUE), 50)
      sample_matrix <- sample_matrix[top_samples, ]
    }

    main_heatmap <- Heatmap(
      sample_matrix,
      name = "MP\nFraction",
      col = fraction_colors,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_names = TRUE,
      show_row_names = FALSE, # Too many sample names
      column_names_side = "bottom",
      column_names_gp = gpar(fontsize = 11, fontface = "bold"),
      top_annotation = HeatmapAnnotation(
        "# Samples" = n_samples_anno,
        "Mean Fraction" = mean_fraction_anno,
        annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
      ),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        grid_height = unit(0.4, "cm")
      )
    )
  }

  # -------------------------------------------------------------------------
  # Panel 4: Bottom Statistics
  # -------------------------------------------------------------------------

  # Total cells per MP (from patient data if available)
  if (!is.null(patient_data)) {
    total_cells <- sapply(mp_names, function(mp) {
      count_col <- paste0(mp, "_count")
      if (count_col %in% colnames(patient_data)) {
        sum(patient_data[[count_col]], na.rm = TRUE)
      } else {
        0
      }
    })
  } else {
    total_cells <- rep(0, length(mp_names))
    names(total_cells) <- mp_names
  }

  cell_count_anno <- anno_barplot(
    total_cells,
    bar_width = 0.8,
    gp = gpar(fill = "#FD8D3C", col = "#D94801"),
    height = unit(1.2, "cm"),
    axis_param = list(gp = gpar(fontsize = 8))
  )

  # -------------------------------------------------------------------------
  # Combine All Panels
  # -------------------------------------------------------------------------

  # Build combined heatmap
  if (!is.null(clinical_heatmap)) {
    combined_heatmap <- clinical_heatmap %v% main_heatmap
  } else {
    combined_heatmap <- main_heatmap
  }

  # Add bottom annotation
  combined_heatmap <- combined_heatmap %v%
    Heatmap(
      matrix(rep(1, length(mp_names)), nrow = 1),
      show_heatmap_legend = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      top_annotation = HeatmapAnnotation(
        "Total Cells" = cell_count_anno,
        annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
      ),
      height = unit(0.1, "cm")
    )

  # =========================================================================
  # PART 6: Save Heatmap
  # =========================================================================

  # Calculate dimensions
  n_mps <- length(mp_names)
  fig_width <- max(10, n_mps * 1.0 + 4)

  if (!is.null(expression_heatmap_data)) {
    fig_height <- max(12, nrow(expression_heatmap_data) * 0.15 + 6)
  } else {
    fig_height <- max(10, min(nrow(sample_matrix), 50) * 0.2 + 6)
  }

  # Create title
  n_samples <- nrow(sample_mp_fractions)
  n_cells <- nrow(clinical_data)
  title_text <- sprintf(
    "Meta-Program Characteristics\n%d Programs | %d Samples | %d Cells",
    length(mp_names), n_samples, n_cells
  )

  # Save PDF
  pdf(file.path(save_path, "mp_features_heatmap_comprehensive.pdf"),
    width = fig_width, height = fig_height
  )

  draw(combined_heatmap,
    column_title = title_text,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    padding = unit(c(3, 3, 3, 3), "mm")
  )

  dev.off()

  cat(sprintf("✅ Comprehensive MP Features Heatmap saved to %s\n", save_path))
  cat("📁 Files: mp_features_heatmap_comprehensive.pdf & .png\n")

  return(NULL)
}

# Helper functions for the proper data structure
calculate_mp_weighted_expression <- function(expression_subset, sample_mp_fractions,
                                             clinical_data, mp_names, genes) {
  # Create MP expression profiles weighted by sample fractions
  mp_expression <- matrix(0, nrow = length(genes), ncol = length(mp_names))
  rownames(mp_expression) <- genes
  colnames(mp_expression) <- mp_names

  # Match samples between expression and fraction data
  sample_ids_expr <- unique(clinical_data$Sample_ID)
  sample_ids_fractions <- sample_mp_fractions$Sample_ID
  common_samples <- intersect(sample_ids_expr, sample_ids_fractions)

  cat(sprintf("     Found %d common samples for expression analysis\n", length(common_samples)))

  for (sample_id in common_samples) {
    # Get cells from this sample
    sample_cells <- clinical_data$V1[clinical_data$Sample_ID == sample_id]
    sample_cells_in_expr <- intersect(sample_cells, colnames(expression_subset))

    if (length(sample_cells_in_expr) > 0) {
      # Get sample expression (mean across cells)
      if (length(sample_cells_in_expr) == 1) {
        sample_expr <- expression_subset[genes, sample_cells_in_expr]
      } else {
        sample_expr <- rowMeans(expression_subset[genes, sample_cells_in_expr], na.rm = TRUE)
      }

      # Get MP fractions for this sample
      sample_fractions <- sample_mp_fractions[sample_mp_fractions$Sample_ID == sample_id, ]

      if (nrow(sample_fractions) > 0) {
        # Weight expression by MP fractions
        for (mp in mp_names) {
          weight <- sample_fractions$MP_fraction[sample_fractions$MP == mp]
          mp_expression[, mp] <- mp_expression[, mp] + sample_expr * weight
        }
      }
    }
  }

  # Normalize by total weights (sample counts)
  for (mp in mp_names) {
    total_weight <- sum(sample_mp_fractions$MP == mp, na.rm = TRUE)
    if (total_weight > 0) {
      mp_expression[, mp] <- mp_expression[, mp] / total_weight
    }
  }

  mp_expression <- apply(mp_expression, 2, function(x) {
    # Z-score normalize each MP
    if (sd(x, na.rm = TRUE) > 0) {
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    } else {
      x
    }
  })

  return(mp_expression)
}

calculate_clinical_associations_from_samples <- function(sample_mp_fractions, mp_names) {
  # Look for clinical variables in sample data
  clinical_vars <- c("Response", "Treatment_Stage", "Microsatellite_Status")
  available_vars <- intersect(clinical_vars, colnames(sample_mp_fractions))

  if (length(available_vars) == 0) {
    cat("   No clinical variables found in sample data\n")
    return(NULL)
  }

  associations <- matrix(0, nrow = length(available_vars), ncol = length(mp_names))
  rownames(associations) <- available_vars
  colnames(associations) <- mp_names

  for (var in available_vars) {
    cat(sprintf("   Calculating sample-level associations for %s...\n", var))

    for (mp in mp_names) {
      rows_ <- sample_mp_fractions$MP == mp
      mp_fractions <- sample_mp_fractions$MP_fraction[rows_]

      clinical_values <- sample_mp_fractions[rows_, ]
      clinical_values <- as.data.frame(clinical_values)[, var]

      # Remove NA and empty values
      valid_idx <- !is.na(mp_fractions) & !is.na(clinical_values) & clinical_values != ""

      if (sum(valid_idx) > 5) {
        mp_clean <- mp_fractions[valid_idx]
        clinical_clean <- clinical_values[valid_idx]

        if (is.character(clinical_clean) || is.factor(clinical_clean)) {
          # Binary categorical variable
          unique_vals <- unique(clinical_clean)
          if (length(unique_vals) == 2) {
            val1 <- unique_vals[1]
            val2 <- unique_vals[2]

            mean1 <- mean(mp_clean[clinical_clean == val1], na.rm = TRUE)
            mean2 <- mean(mp_clean[clinical_clean == val2], na.rm = TRUE)

            # Effect size (Cohen's d)
            pooled_sd <- sd(mp_clean, na.rm = TRUE)
            if (pooled_sd > 0) {
              associations[var, mp] <- (mean1 - mean2) / pooled_sd
            }
          }
        } else {
          # Numeric variable
          clinical_numeric <- as.numeric(clinical_clean)
          if (sd(clinical_numeric, na.rm = TRUE) > 0) {
            associations[var, mp] <- cor(mp_clean, clinical_numeric, use = "complete.obs")
          }
        }
      }
    }
  }

  return(as.data.frame(associations))
}

calculate_sample_level_stats <- function(sample_mp_fractions, mp_names) {
  stats <- list()

  # Number of samples with non-zero fraction for each MP
  stats$n_samples <- sapply(mp_names, function(mp) {
    sum(sample_mp_fractions$MP_fraction[sample_mp_fractions$MP == mp] > 0.01, na.rm = TRUE) # >1% threshold
  })

  # Mean fraction across all samples
  stats$mean_fraction <- sapply(mp_names, function(mp) {
    mean(sample_mp_fractions$MP_fraction[sample_mp_fractions$MP == mp], na.rm = TRUE)
  })

  # Max fraction observed
  stats$max_fraction <- sapply(mp_names, function(mp) {
    max(sample_mp_fractions$MP_fraction[sample_mp_fractions$MP == mp], na.rm = TRUE)
  })

  names(stats$n_samples) <- mp_names
  names(stats$mean_fraction) <- mp_names
  names(stats$max_fraction) <- mp_names

  return(stats)
}

create_gene_annotation_from_markers <- function(genes, top_genes_list, mp_names) {
  gene_mp_map <- rep("Unassigned", length(genes))
  names(gene_mp_map) <- genes

  # Assign each gene to its primary MP
  for (mp in mp_names) {
    if (mp %in% names(top_genes_list)) {
      mp_genes <- top_genes_list[[mp]]
      for (gene in mp_genes) {
        if (gene %in% genes && gene_mp_map[gene] == "Unassigned") {
          gene_mp_map[gene] <- mp
        }
      }
    }
  }

  return(gene_mp_map)
}

create_mp_color_scheme <- function(mp_names) {
  n_mps <- length(mp_names)

  if (n_mps <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n_mps), "Set2")[1:n_mps]
  } else if (n_mps <= 12) {
    colors <- RColorBrewer::brewer.pal(n_mps, "Set3")
  } else {
    colors <- rainbow(n_mps, end = 0.85)
  }

  names(colors) <- mp_names
  colors["Unassigned"] <- "#CCCCCC"

  return(colors)
}
