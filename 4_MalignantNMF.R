# MSI-H CRC Meta-Program Resistance Analysis - R Script

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(viridis)
  library(RColorBrewer)
  library(corrplot)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(scales)
  library(readr)
  library(dplyr)
  library(reshape2)
  library(broom)
  library(ggrepel)
  library(cowplot)
})

# Set working directory (adjust path as needed)
workDir <- "/mnt/public/lyx/CRC_dMMR"
save_path <- file.path(workDir, "Results", "integration_results")
result_path <- file.path(workDir, "Results", "NMF_results")
figurePath <- file.path(workDir, "Figures", "NMF_results")

source("malignant_functions.r")

# Define consistent color palettes
colors <- list(
  # Response colors
  response = c("pCR" = "#2E8B57", "non_pCR" = "#DC143C"), # SeaGreen vs Crimson

  # Treatment stage colors
  stage = c("Pre" = "#4682B4", "Post" = "#FF6347"), # SteelBlue vs Tomato

  # Microsatellite status colors
  ms_status = c("MSI" = "#FF8C00", "MSS" = "#9370DB"), # DarkOrange vs MediumPurple

  # Treatment strategy colors
  treatment = c(
    "Anti-PD1" = "#32CD32",
    "Anti-PD1 plus Celecoxib" = "#20B2AA",
    "Anti-PD1 plus CapeOx" = "#4169E1"
  ),

  # MP classification colors
  mp_class = c(
    "MSI_enriched" = "#FF4500", # OrangeRed
    "MSI_slight" = "#FFA500", # Orange
    "Shared" = "#808080", # Gray
    "MSS_slight" = "#9932CC", # DarkOrchid
    "MSS_enriched" = "#4B0082" # Indigo
  ),

  # Significance colors
  significance = c("Significant" = "#B22222", "Non-significant" = "#696969")
)

cat("🧬 MSI-H CRC Meta-Program Resistance Analysis - R Script\n")

# =============================================================================
# SECTION 1: DATA LOADING AND PREPROCESSING
# =============================================================================

cat("\n📂 Loading Data Files...\n")

# Load all sample MP fraction data
all_sample_data <- read_csv(file.path(result_path, "all_sample_mp_fraction_data.csv"), show_col_types = FALSE)
cat("✅ Loaded all sample MP fraction data:", nrow(all_sample_data), "rows\n")

# Load task-specific data
task1_data <- file.path(result_path, "task1_intrinsic_resistance_sample_data.csv")
task2_data <- file.path(result_path, "task2_acquired_resistance_sample_data.csv")
task3_data <- file.path(result_path, "task3_msi_mss_similarity_sample_data.csv")

if (file.exists(task1_data)) {
  task1_data <- read_csv(task1_data, show_col_types = FALSE)
  cat("✅ Loaded Task 1 data:", nrow(task1_data), "rows\n")
}

if (file.exists(task2_data)) {
  task2_data <- read_csv(task2_data, show_col_types = FALSE)
  cat("✅ Loaded Task 2 data:", nrow(task2_data), "rows\n")
}

if (file.exists(task3_data)) {
  task3_data <- read_csv(task3_data, show_col_types = FALSE)
  cat("✅ Loaded Task 3 data:", nrow(task3_data), "rows\n")
}

# Load statistical results
task1_stats <- file.path(result_path, "task1_intrinsic_resistance_statistical_results.csv")
task2_stats <- file.path(result_path, "task2_acquired_resistance_statistical_results.csv")
task3_stats <- file.path(result_path, "task3_msi_mss_similarity_statistical_results.csv")

if (file.exists(task1_stats)) {
  task1_stats <- read_csv(task1_stats, show_col_types = FALSE)
}

if (file.exists(task2_stats)) {
  task2_stats <- read_csv(task2_stats, show_col_types = FALSE)
}

if (file.exists(task3_stats)) {
  task3_stats <- read_csv(task3_stats, show_col_types = FALSE)
}

# Load MP classification if available
mp_classification <- file.path(result_path, "task3_mp_classification.csv")
if (file.exists(mp_classification)) {
  mp_classification <- read_csv(mp_classification, show_col_types = FALSE)
}

# =============================================================================
# SECTION 2: DATA OVERVIEW AND QUALITY CHECKS
# =============================================================================

cat("\n📊 Data Overview and Quality Checks...\n")

# Create data overview plot
data_overview <- all_sample_data %>%
  group_by(Microsatellite_Status, Treatment_Stage, Response) %>%
  summarise(
    n_samples = n_distinct(Sample_ID),
    n_patients = n_distinct(Patient_ID),
    .groups = "drop"
  ) %>%
  filter(!is.na(Microsatellite_Status), !is.na(Response))

overview_plot <- data_overview %>%
  ggplot(aes(x = Treatment_Stage, y = n_samples, fill = Response)) +
  geom_col(position = "dodge", alpha = 0.8) +
  facet_wrap(~Microsatellite_Status, scales = "free_y") +
  scale_fill_manual(values = colors$response) +
  labs(
    title = "Data Overview: Sample Distribution",
    subtitle = "Number of samples by microsatellite status, treatment stage, and response",
    x = "Treatment Stage",
    y = "Number of Samples",
    fill = "Response"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

save_plot(overview_plot, figurePath, "01_data_overview", width = 12, height = 8)

# =============================================================================
# SECTION 3: TASK 1 - INTRINSIC RESISTANCE ANALYSIS
# =============================================================================

if (!is.null(task1_data)) {
  cat("\n🎯 Task 1: Intrinsic Resistance Analysis (MSI Pre-treatment)\n")

  # 1. Boxplot comparison: pCR vs non_pCR
  task1_boxplot <- task1_data %>%
    ggplot(aes(x = Response, y = MP_fraction, fill = Response)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    facet_wrap(~MP, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = colors$response) +
    labs(
      title = "Task 1: Intrinsic Resistance - MP Fractions by Response",
      subtitle = "MSI tumors, pre-treatment samples only",
      x = "Response",
      y = "MP Fraction",
      fill = "Response"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    stat_compare_means(
      aes(label = ..p..),
      method = "t.test",
      comparisons = list(c("pCR", "non_pCR")),
      size = 5
    )

  save_plot(task1_boxplot, figurePath, "02_task1_intrinsic_resistance_boxplot", width = 16, height = 12)

  # 2. Statistical significance volcano plot
  if (!is.null(task1_stats)) {
    task1_volcano <- task1_stats %>%
      mutate(
        significance = add_significance(p_corrected),
        is_significant = p_corrected < 0.05 & abs(cohens_d) > 0.5,
        direction = case_when(
          cohens_d > 0.5 & p_corrected < 0.05 ~ "Higher in non_pCR",
          cohens_d < -0.5 & p_corrected < 0.05 ~ "Higher in pCR",
          TRUE ~ "Non-significant"
        )
      ) %>%
      ggplot(aes(x = cohens_d, y = -log10(p_corrected))) +
      geom_point(aes(color = direction, size = abs(mean_difference)), alpha = 0.7) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.6) +
      geom_text_repel(
        data = . %>% filter(is_significant),
        aes(label = MP),
        size = 3,
        max.overlaps = 20
      ) +
      scale_color_manual(
        values = c(
          "Higher in non_pCR" = colors$response[["non_pCR"]],
          "Higher in pCR" = colors$response[["pCR"]],
          "Non-significant" = "gray50"
        )
      ) +
      scale_size_continuous(range = c(1, 4)) +
      labs(
        title = "Task 1: Intrinsic Resistance - Statistical Significance",
        subtitle = "Effect size vs significance (FDR-corrected p-value)",
        x = "Cohen's d (Effect Size)",
        y = "-log10(FDR-corrected p-value)",
        color = "Direction",
        size = "Mean Difference"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom"
      )

    save_plot(task1_volcano, figurePath, "03_task1_intrinsic_volcano_plot", width = 12, height = 10)
  }

  # 3. Effect size heatmap
  if (!is.null(task1_stats)) {
    task1_heatmap_data <- task1_stats %>%
      select(MP, cohens_d, p_corrected) %>%
      mutate(
        significance = add_significance(p_corrected),
        cohens_d_capped = pmax(-2, pmin(2, cohens_d)) # Cap for visualization
      )

    task1_heatmap <- task1_heatmap_data %>%
      ggplot(aes(x = 1, y = reorder(MP, cohens_d), fill = cohens_d_capped)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = p_corrected), color = "black", size = 3, fontface = "bold") +
      scale_fill_gradient2(
        low = colors$response[["pCR"]],
        mid = "white",
        high = colors$response[["non_pCR"]],
        midpoint = 0,
        limits = c(-2, 2),
        name = "Cohen's d"
      ) +
      labs(
        title = "Task 1: Intrinsic Resistance - Effect Size Heatmap",
        subtitle = "Positive values: higher in non_pCR; Negative values: higher in pCR",
        x = "",
        y = "Meta-Program"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )

    save_plot(task1_heatmap, figurePath, "04_task1_intrinsic_effect_size_heatmap", width = 6, height = 8)
  }
}

# =============================================================================
# SECTION 4: TASK 2 - ACQUIRED RESISTANCE ANALYSIS
# =============================================================================

if (!is.null(task2_data)) {
  cat("\n🔥 Task 2: Acquired Resistance Analysis (MSI non_pCR Pre vs Post)\n")

  # 1. Paired comparison plot
  task2_paired_plot <- task2_data %>%
    ggplot(aes(x = Treatment_Stage, y = MP_fraction, color = Treatment_Stage)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_line(aes(group = Patient_ID), alpha = 0.4, size = 0.5) +
    geom_point(alpha = 0.6, size = 1.5) +
    facet_wrap(~MP, scales = "free_y", ncol = 4) +
    scale_color_manual(values = colors$stage) +
    labs(
      title = "Task 2: Acquired Resistance - MP Changes During Treatment",
      subtitle = "MSI non_pCR patients, paired Pre vs Post samples",
      x = "Treatment Stage",
      y = "MP Fraction",
      color = "Treatment Stage"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    stat_compare_means(
      method = "t.test",
      size = 5
    )

  save_plot(task2_paired_plot, figurePath, "05_task2_acquired_resistance_paired", width = 16, height = 12)

  # 2. Change direction plot (Post - Pre)
  task2_changes <- task2_data %>%
    select(Patient_ID, MP, MP_fraction, Treatment_Stage) %>%
    pivot_wider(names_from = Treatment_Stage, values_from = MP_fraction) %>%
    filter(!is.na(Pre), !is.na(Post)) %>%
    mutate(change = Post - Pre)

  if (nrow(task2_changes) > 0) {
    task2_change_plot <- task2_changes %>%
      ggplot(aes(x = MP, y = change)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
      geom_boxplot(aes(fill = ifelse(change > 0, "Increased", "Decreased")), alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      scale_fill_manual(
        values = c("Increased" = colors$stage[["Post"]], "Decreased" = colors$stage[["Pre"]])
      ) +
      labs(
        title = "Task 2: Treatment-Induced MP Changes",
        subtitle = "Change in MP fraction (Post - Pre) in MSI non_pCR patients",
        x = "Meta-Program",
        y = "Change in MP Fraction (Post - Pre)",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )

    save_plot(task2_change_plot, figurePath, "06_task2_treatment_induced_changes", width = 12, height = 8)
  }

  # 3. Statistical significance for acquired resistance
  if (!is.null(task2_stats)) {
    task2_volcano <- task2_stats %>%
      mutate(
        significance = add_significance(p_corrected),
        is_significant = p_corrected < 0.05 & abs(cohens_d) > 0.5,
        direction = case_when(
          cohens_d > 0.5 & p_corrected < 0.05 ~ "Increased Post-treatment",
          cohens_d < -0.5 & p_corrected < 0.05 ~ "Decreased Post-treatment",
          TRUE ~ "Non-significant"
        )
      ) %>%
      ggplot(aes(x = cohens_d, y = -log10(p_corrected))) +
      geom_point(aes(color = direction, size = abs(mean_difference)), alpha = 0.7) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.6) +
      geom_text_repel(
        data = . %>% filter(is_significant),
        aes(label = MP),
        size = 3,
        max.overlaps = 20
      ) +
      scale_color_manual(
        values = c(
          "Increased Post-treatment" = colors$stage[["Post"]],
          "Decreased Post-treatment" = colors$stage[["Pre"]],
          "Non-significant" = "gray50"
        )
      ) +
      scale_size_continuous(range = c(1, 4)) +
      labs(
        title = "Task 2: Acquired Resistance - Statistical Significance",
        subtitle = "Effect size vs significance for treatment-induced changes",
        x = "Cohen's d (Effect Size)",
        y = "-log10(FDR-corrected p-value)",
        color = "Direction",
        size = "Mean Difference"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom"
      )

    save_plot(task2_volcano, figurePath, "07_task2_acquired_volcano_plot", width = 12, height = 10)
  }
}

# =============================================================================
# SECTION 5: TASK 3 - MSI VS MSS SIMILARITY ANALYSIS
# =============================================================================

if (!is.null(task3_data)) {
  cat("\n⚡ Task 3: MSI vs MSS Similarity Analysis\n")

  # 1. MP fraction comparison: MSI vs MSS
  task3_comparison <- task3_data %>%
    ggplot(aes(x = Microsatellite_Status, y = MP_fraction, fill = Microsatellite_Status)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    facet_wrap(~MP, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = colors$ms_status) +
    labs(
      title = "Task 3: MSI vs MSS Meta-Program Comparison",
      subtitle = "MP fractions in MSI vs MSS treated tumors",
      x = "Microsatellite Status",
      y = "MP Fraction",
      fill = "Microsatellite Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    stat_compare_means(
      aes(label = "p"),
      method = "t.test",
      comparisons = list(c("MSI", "MSS")),
      size = 5
    )

  save_plot(task3_comparison, figurePath, "08_task3_msi_mss_comparison", width = 16, height = 12)

  # 2. Correlation scatter plot
  task3_correlation_data <- task3_data %>%
    group_by(MP, Microsatellite_Status) %>%
    summarise(mean_fraction = mean(MP_fraction), .groups = "drop") %>%
    pivot_wider(names_from = Microsatellite_Status, values_from = mean_fraction) %>%
    filter(!is.na(MSI), !is.na(MSS))

  if (nrow(task3_correlation_data) > 0) {
    correlation_coef <- cor(task3_correlation_data$MSI, task3_correlation_data$MSS,
      use = "complete.obs"
    )

    task3_scatter <- task3_correlation_data %>%
      ggplot(aes(x = MSI, y = MSS)) +
      geom_point(size = 3, alpha = 0.7, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
      geom_text_repel(aes(label = MP), size = 3, max.overlaps = 15) +
      labs(
        title = "Task 3: MSI vs MSS Meta-Program Correlation",
        subtitle = paste0("Correlation coefficient: r = ", round(correlation_coef, 3)),
        x = "Mean MP Fraction in MSI",
        y = "Mean MP Fraction in MSS"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )

    save_plot(task3_scatter, figurePath, "09_task3_msi_mss_correlation", width = 10, height = 8)
  }

  # 3. MP classification plot
  if (!is.null(mp_classification)) {
    mp_class_plot <- mp_classification %>%
      filter(!is.na(category)) %>%
      mutate(
        category = factor(category, levels = c(
          "MSI_enriched", "MSI_slight", "Shared",
          "MSS_slight", "MSS_enriched"
        )),
        is_significant = significant & abs(cohens_d) > 0.5
      ) %>%
      ggplot(aes(x = cohens_d, y = reorder(MP, cohens_d), color = category)) +
      geom_vline(xintercept = c(-0.5, -0.2, 0.2, 0.5), linetype = "dashed", alpha = 0.4) +
      geom_point(aes(size = is_significant, alpha = is_significant)) +
      scale_color_manual(values = colors$mp_class) +
      scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.5)) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6)) +
      labs(
        title = "Task 3: Meta-Program Classification",
        subtitle = "MSI-enriched vs MSS-enriched vs Shared MPs",
        x = "Cohen's d (MSS - MSI)",
        y = "Meta-Program",
        color = "Classification",
        size = "Significant",
        alpha = "Significant"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom"
      )

    save_plot(mp_class_plot, figurePath, "10_task3_mp_classification", width = 12, height = 10)
  }
}

# =============================================================================
# SECTION 6: COMPREHENSIVE SUMMARY PLOTS
# =============================================================================

cat("\n📊 Creating Comprehensive Summary Plots...\n")

# 1. Combined effect size comparison across all tasks
combined_effects <- data.frame()

if (!is.null(task1_stats)) {
  task1_effects <- task1_stats %>%
    mutate(
      Task = "Task 1: Intrinsic Resistance",
      Comparison = "non_pCR vs pCR"
    ) %>%
    select(Task, Comparison, MP, cohens_d, p_corrected)
  combined_effects <- rbind(combined_effects, task1_effects)
}

if (!is.null(task2_stats)) {
  task2_effects <- task2_stats %>%
    mutate(
      Task = "Task 2: Acquired Resistance",
      Comparison = "Post vs Pre"
    ) %>%
    select(Task, Comparison, MP, cohens_d, p_corrected)
  combined_effects <- rbind(combined_effects, task2_effects)
}

if (!is.null(task3_stats)) {
  task3_effects <- task3_stats %>%
    mutate(
      Task = "Task 3: MSI vs MSS",
      Comparison = "MSS vs MSI"
    ) %>%
    select(Task, Comparison, MP, cohens_d, p_corrected)
  combined_effects <- rbind(combined_effects, task3_effects)
}

if (nrow(combined_effects) > 0) {
  combined_effects_plot <- combined_effects %>%
    mutate(
      significance = add_significance(p_corrected),
      is_significant = p_corrected < 0.05 & abs(cohens_d) > 0.5
    ) %>%
    ggplot(aes(x = MP, y = cohens_d, fill = Task)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(
      aes(label = ifelse(is_significant, significance, "")),
      position = position_dodge(width = 0.9),
      vjust = ifelse(combined_effects$cohens_d > 0, -0.3, 1.3),
      size = 3
    ) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.6) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = "Comprehensive Effect Size Comparison Across All Tasks",
      subtitle = "Cohen's d values for different resistance analyses",
      x = "Meta-Program",
      y = "Cohen's d (Effect Size)",
      fill = "Analysis Task"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  save_plot(combined_effects_plot, figurePath, "11_comprehensive_effect_sizes", width = 16, height = 10)
}

# 2. Summary statistics table plot
if (nrow(combined_effects) > 0) {
  summary_stats <- combined_effects %>%
    group_by(Task) %>%
    summarise(
      Total_MPs = n(),
      Significant_MPs = sum(p_corrected < 0.05 & abs(cohens_d) > 0.5),
      Large_Effects = sum(abs(cohens_d) > 0.8),
      Mean_Effect_Size = round(mean(abs(cohens_d)), 3),
      .groups = "drop"
    )

  # Convert to long format for plotting
  summary_long <- summary_stats %>%
    pivot_longer(cols = -Task, names_to = "Metric", values_to = "Value")

  summary_table_plot <- summary_long %>%
    ggplot(aes(x = Task, y = Metric, fill = Value)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = Value), color = "white", size = 4, fontface = "bold") +
    scale_fill_viridis_c(name = "Count/Value") +
    labs(
      title = "Summary Statistics by Analysis Task",
      x = "Analysis Task",
      y = "Metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  save_plot(summary_table_plot, figurePath, "12_summary_statistics_table", width = 12, height = 6)
}

# =============================================================================
# SECTION 7: EXPORT SUMMARY STATISTICS
# =============================================================================

cat("\n💾 Exporting Summary Statistics...\n")

# Create comprehensive summary report
if (exists("summary_stats")) {
  write_csv(summary_stats, file.path(result_path, "summary_statistics.csv"))
  cat("✅ Exported summary statistics\n")
}

# Export statistical results with R calculations
if (nrow(combined_effects) > 0) {
  r_statistical_summary <- combined_effects %>%
    group_by(Task, MP) %>%
    summarise(
      cohen_d = first(cohens_d),
      p_value = first(p_corrected),
      significance = add_significance(first(p_corrected)),
      effect_magnitude = case_when(
        abs(first(cohens_d)) >= 0.8 ~ "Large",
        abs(first(cohens_d)) >= 0.5 ~ "Medium",
        abs(first(cohens_d)) >= 0.2 ~ "Small",
        TRUE ~ "Negligible"
      ),
      is_significant = first(p_corrected) < 0.05 & abs(first(cohens_d)) > 0.5,
      .groups = "drop"
    )

  write_csv(r_statistical_summary, file.path(result_path, "r_statistical_summary.csv"))
  cat("✅ Exported R statistical summary\n")
}

# =============================================================================
# SECTION 8: ADVANCED VISUALIZATIONS
# =============================================================================

cat("\n🎨 Creating Advanced Visualizations...\n")

# 1. MP Fraction Heatmap across all conditions
if (!is.null(all_sample_data)) {
  # Create sample-level summary for heatmap
  heatmap_data <- all_sample_data %>%
    filter(!is.na(Microsatellite_Status), !is.na(Response)) %>%
    group_by(MP, Microsatellite_Status, Response, Treatment_Stage) %>%
    summarise(mean_fraction = mean(MP_fraction), .groups = "drop") %>%
    unite("Condition", Microsatellite_Status, Response, Treatment_Stage, sep = "_") %>%
    pivot_wider(names_from = Condition, values_from = mean_fraction, values_fill = 0)

  heatmap_data <- heatmap_data[, !endsWith(colnames(heatmap_data), suffix = "On")]

  if (nrow(heatmap_data) > 1 && ncol(heatmap_data) > 2) {
    # Prepare matrix for heatmap
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    rownames(heatmap_matrix) <- heatmap_data$MP

    # Create annotation for conditions
    condition_annotation <- data.frame()
    anno_ <- strsplit(colnames(heatmap_matrix), "_")
    for (i in seq_along(anno_)) {
      if (length(anno_[[i]]) == 3) {
        temp_ <- setNames(anno_[[i]], c("MS_Status", "Response", "Stage"))
      }
      if (length(anno_[[i]]) == 4) {
        temp_ <- setNames(c(anno_[[i]][1], paste0(anno_[[i]][2], "_", anno_[[i]][3]), anno_[[i]][4]), c("MS_Status", "Response", "Stage"))
      }
      condition_annotation <- rbind(condition_annotation, temp_)
    }
    colnames(condition_annotation) <- c("MS_Status", "Response", "Stage")
    rownames(condition_annotation) <- colnames(heatmap_matrix)


    # Define annotation colors
    annotation_colors <- list(
      MS_Status = colors$ms_status[condition_annotation$MS_Status],
      Response = colors$response[condition_annotation$Response],
      Stage = colors$stage[condition_annotation$Stage]
    )

    # Create heatmap
    pdf(file.path(figurePath, "13_comprehensive_mp_heatmap.pdf"), width = 12, height = 10)
    pheatmap(
      heatmap_matrix,
      annotation_col = condition_annotation,
      annotation_colors = annotation_colors,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      scale = "row",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Meta-Program Fractions Across All Conditions",
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8
    )
    dev.off()

    cat("✅ Saved comprehensive MP heatmap\n")
  }
}

# 2. Multi-panel figure combining key results
if (!is.null(task1_data) && !is.null(task3_data)) {
  # Panel A: Task 1 key result
  panel_a <- task1_data %>%
    group_by(MP, Response) %>%
    summarise(mean_fraction = mean(MP_fraction), .groups = "drop") %>%
    ggplot(aes(x = reorder(MP, mean_fraction), y = mean_fraction, fill = Response)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = colors$response) +
    labs(
      title = "A) Intrinsic Resistance",
      x = "Meta-Program",
      y = "Mean MP Fraction",
      fill = "Response"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  # Panel B: Task 3 key result
  panel_b <- task3_data %>%
    group_by(MP, Microsatellite_Status) %>%
    summarise(mean_fraction = mean(MP_fraction), .groups = "drop") %>%
    ggplot(aes(x = reorder(MP, mean_fraction), y = mean_fraction, fill = Microsatellite_Status)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = colors$ms_status) +
    labs(
      title = "B) MSI vs MSS Patterns",
      x = "Meta-Program",
      y = "Mean MP Fraction",
      fill = "Microsatellite Status"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  # Combine panels
  multi_panel <- panel_a / panel_b +
    plot_annotation(
      title = "Key Findings: MSI-H CRC Meta-Program Resistance Analysis",
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
    )

  save_plot(multi_panel, figurePath, "14_multi_panel_key_findings", width = 14, height = 12)
}

# 3. Treatment strategy specific analysis
if (!is.null(all_sample_data)) {
  treatment_analysis <- all_sample_data %>%
    filter(
      !is.na(Treatment_Strategy),
      Treatment_Strategy %in% c("Anti-PD1", "Anti-PD1 plus Celecoxib", "Anti-PD1 plus CapeOx"),
      Microsatellite_Status == "MSI"
    ) %>%
    group_by(MP, Treatment_Strategy, Response) %>%
    summarise(
      mean_fraction = mean(MP_fraction),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 3) # Only include groups with sufficient samples

  if (nrow(treatment_analysis) > 0) {
    treatment_plot <- treatment_analysis %>%
      ggplot(aes(x = MP, y = mean_fraction, fill = interaction(Treatment_Strategy, Response))) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(
        values = c(
          "Anti-PD1.pCR" = "#90EE90",
          "Anti-PD1.non_pCR" = "#FF6B6B",
          "Anti-PD1 plus Celecoxib.pCR" = "#87CEEB",
          "Anti-PD1 plus Celecoxib.non_pCR" = "#DDA0DD",
          "Anti-PD1 plus CapeOx.pCR" = "#98FB98",
          "Anti-PD1 plus CapeOx.non_pCR" = "#F0E68C"
        ),
        name = "Treatment & Response"
      ) +
      labs(
        title = "Treatment Strategy Analysis",
        subtitle = "MP fractions by treatment strategy and response (MSI tumors)",
        x = "Meta-Program",
        y = "Mean MP Fraction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      ) +
      guides(fill = guide_legend(ncol = 3))

    save_plot(treatment_plot, figurePath, "15_treatment_strategy_analysis", width = 16, height = 10)
  }
}

# =============================================================================
# SECTION 9: STATISTICAL VALIDATION AND REPRODUCIBILITY
# =============================================================================

cat("\n🔬 Statistical Validation and Reproducibility Checks...\n")

# Reproduce key statistical tests from Python results
reproducibility_results <- list()

# Task 1 validation
if (!is.null(task1_data) && !is.null(task1_stats)) {
  cat("🔍 Validating Task 1 statistics...\n")

  task1_validation <- task1_data %>%
    group_by(MP) %>%
    do({
      data <- .
      pcr_vals <- data$MP_fraction[data$Response == "pCR"]
      nonpcr_vals <- data$MP_fraction[data$Response == "non_pCR"]

      if (length(pcr_vals) >= 3 && length(nonpcr_vals) >= 3) {
        test_result <- t.test(nonpcr_vals, pcr_vals)

        # Calculate Cohen's d
        pooled_sd <- sqrt(((length(nonpcr_vals) - 1) * var(nonpcr_vals) +
          (length(pcr_vals) - 1) * var(pcr_vals)) /
          (length(nonpcr_vals) + length(pcr_vals) - 2))
        cohens_d <- (mean(nonpcr_vals) - mean(pcr_vals)) / pooled_sd

        data.frame(
          MP = unique(data$MP),
          r_t_statistic = test_result$statistic,
          r_p_value = test_result$p.value,
          r_cohens_d = cohens_d,
          r_mean_difference = mean(nonpcr_vals) - mean(pcr_vals)
        )
      } else {
        data.frame()
      }
    }) %>%
    ungroup()

  # Compare with Python results
  if (nrow(task1_validation) > 0) {
    task1_comparison <- task1_stats %>%
      select(MP, t_statistic, p_value, cohens_d, mean_difference) %>%
      inner_join(task1_validation, by = "MP") %>%
      mutate(
        t_stat_diff = abs(t_statistic - r_t_statistic),
        cohens_d_diff = abs(cohens_d - r_cohens_d),
        validation_status = ifelse(t_stat_diff < 0.01 & cohens_d_diff < 0.01, "✅ Validated", "⚠️ Discrepancy")
      )

    reproducibility_results$task1 <- task1_comparison
    write_csv(task1_comparison, file.path(result_path, "task1_validation_comparison.csv"))

    cat(
      "✅ Task 1 validation completed -",
      sum(task1_comparison$validation_status == "✅ Validated"), "out of",
      nrow(task1_comparison), "tests validated\n"
    )
  }
}

# Create validation summary plot
if (length(reproducibility_results) > 0 && !is.null(reproducibility_results$task1)) {
  validation_plot <- reproducibility_results$task1 %>%
    ggplot(aes(x = cohens_d, y = r_cohens_d, color = validation_status)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("✅ Validated" = "green", "⚠️ Discrepancy" = "red")) +
    labs(
      title = "Statistical Validation: Python vs R Results",
      subtitle = "Comparison of Cohen's d effect sizes",
      x = "Python Cohen's d",
      y = "R Cohen's d",
      color = "Validation Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom"
    )

  save_plot(validation_plot, figurePath, "16_statistical_validation", width = 10, height = 8)
}

# =============================================================================
# SECTION 10: FINAL SUMMARY AND EXPORT
# =============================================================================

cat("\n📋 Creating Final Summary Report...\n")

# Generate comprehensive analysis report
analysis_summary <- list(
  analysis_date = Sys.Date(),
  total_samples = ifelse(!is.null(all_sample_data), length(unique(all_sample_data$Sample_ID)), 0),
  total_patients = ifelse(!is.null(all_sample_data), length(unique(all_sample_data$Patient_ID)), 0),
  tasks_completed = c(
    task1 = !is.null(task1_data),
    task2 = !is.null(task2_data),
    task3 = !is.null(task3_data)
  ),
  figures_generated = list.files(figurePath, pattern = "\\.(pdf|png)$"),
  key_findings = list()
)

# Add key findings
if (!is.null(task1_stats)) {
  significant_intrinsic <- task1_stats %>%
    filter(p_corrected < 0.05 & abs(cohens_d) > 0.5) %>%
    pull(MP)
  analysis_summary$key_findings$intrinsic_resistance_mps <- significant_intrinsic
}

if (!is.null(task2_stats)) {
  significant_acquired <- task2_stats %>%
    filter(p_corrected < 0.05 & cohens_d > 0.5) %>%
    pull(MP)
  analysis_summary$key_findings$acquired_resistance_mps <- significant_acquired
}

if (!is.null(mp_classification)) {
  msi_enriched <- mp_classification %>%
    filter(category == "MSI_enriched" & significant) %>%
    pull(MP)
  mss_enriched <- mp_classification %>%
    filter(category == "MSS_enriched" & significant) %>%
    pull(MP)

  analysis_summary$key_findings$msi_enriched_mps <- msi_enriched
  analysis_summary$key_findings$mss_enriched_mps <- mss_enriched
}


# # Save analysis summary
# saveRDS(analysis_summary, "R_figures/analysis_summary.rds")

# Create final summary text file
summary_text <- c(
  "MSI-H CRC Meta-Program Resistance Analysis - Final Summary",
  "================================================================",
  "",
  paste("Analysis Date:", Sys.Date()),
  paste("Total Samples Analyzed:", analysis_summary$total_samples),
  paste("Total Patients Analyzed:", analysis_summary$total_patients),
  "",
  "Tasks Completed:",
  paste("  Task 1 (Intrinsic Resistance):", ifelse(analysis_summary$tasks_completed[1], "✅", "❌")),
  paste("  Task 2 (Acquired Resistance):", ifelse(analysis_summary$tasks_completed[2], "✅", "❌")),
  paste("  Task 3 (MSI vs MSS Similarity):", ifelse(analysis_summary$tasks_completed[3], "✅", "❌")),
  "",
  paste("Total Figures Generated:", length(analysis_summary$figures_generated)),
  "",
  "Key Findings:",
  paste("  Intrinsic Resistance MPs:", ifelse(length(analysis_summary$key_findings$intrinsic_resistance_mps) > 0,
    paste(analysis_summary$key_findings$intrinsic_resistance_mps, collapse = ", "),
    "None found"
  )),
  paste("  Acquired Resistance MPs:", ifelse(length(analysis_summary$key_findings$acquired_resistance_mps) > 0,
    paste(analysis_summary$key_findings$acquired_resistance_mps, collapse = ", "),
    "None found"
  )),
  paste("  MSI-Enriched MPs:", ifelse(length(analysis_summary$key_findings$msi_enriched_mps) > 0,
    paste(analysis_summary$key_findings$msi_enriched_mps, collapse = ", "),
    "None found"
  )),
  paste("  MSS-Enriched MPs:", ifelse(length(analysis_summary$key_findings$mss_enriched_mps) > 0,
    paste(analysis_summary$key_findings$mss_enriched_mps, collapse = ", "),
    "None found"
  )),
  "",
  "Files Generated:",
  paste("  -", analysis_summary$figures_generated)
)

writeLines(summary_text, file.path(result_path, "ANALYSIS_SUMMARY.txt"))

cat("\n🎉 ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("📁 All results saved in: ", result_path, "\n")

# Final message with color scheme information
cat("\n🎨 Color Scheme Used (for consistency):\n")
cat("  Response: pCR =", colors$response[["pCR"]], ", non_pCR =", colors$response[["non_pCR"]], "\n")
cat("  Treatment Stage: Pre =", colors$stage[["Pre"]], ", Post =", colors$stage[["Post"]], "\n")
cat("  Microsatellite: MSI =", colors$ms_status[["MSI"]], ", MSS =", colors$ms_status[["MSS"]], "\n")
cat("\n✨ Analysis pipeline completed! Check ", figurePath, "/ directory for all outputs.\n")

## NEW ADD
data_list <- load_sparse_expression_data(file.path(result_path, "NMF_results_for_R"))
# 1. MP Composition Stacked Barplot
if (!is.null(data_list$sample_mp_fractions)) {
  cat("📊 Creating MP composition stacked barplot...\n")
  create_mp_composition_barplot_existing_style(data_list, figurePath)
}

# 2. MP Marker Gene Heatmap
if (!is.null(data_list$mp_marker_genes) && !is.null(data_list$expression_subset)) {
  cat("🔥 Creating MP marker gene heatmap...\n")
  create_mp_composition_barplot_modern(data_list, figurePath)
}

# Comprehensive heatmap
results <- create_mp_features_heatmap(
  data_list = data_list,
  save_path = figurePath,
  top_genes_per_mp = 25
)
