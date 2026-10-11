###############################################################################
# AD-DEG-Analysis.R
# DESeq2 differential expression analysis
# Comparison: Alzheimer vs control
# Covariates: age at death, RIN score, post-mortem interval, and sex
#
# Inputs (place in the working directory, or update the paths):
#   - metadata_filtered.txt  (sample metadata; needs columns: group,
#                             diagnosis, sex, age_at_death, rin_score,
#                             post_mortem_interval)
#   - gene_counts_brain.txt  (raw gene counts, genes x samples)
###############################################################################


# =============================================================================
# 1. Load packages
# =============================================================================

library(dplyr)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(grid)
library(circlize)

# =============================================================================
# 2. Load metadata
# =============================================================================

metadata_filtered <- read.delim(
  "metadata_filtered.txt",
  header = TRUE,
  row.names = 1
)

rownames(metadata_filtered) <- metadata_filtered$group


# =============================================================================
# 3. Prepare metadata and covariates
# =============================================================================

# Convert covariates to numeric.
metadata_filtered$age_at_death <- as.numeric(
  as.character(metadata_filtered$age_at_death)
)
metadata_filtered$post_mortem_interval <- as.numeric(
  as.character(metadata_filtered$post_mortem_interval)
)
metadata_filtered$rin_score <- as.numeric(
  as.character(metadata_filtered$rin_score)
)

# Impute missing covariate values with the corresponding mean.
metadata_filtered$age_at_death[is.na(metadata_filtered$age_at_death)] <-
  mean(metadata_filtered$age_at_death, na.rm = TRUE)

metadata_filtered$post_mortem_interval[is.na(metadata_filtered$post_mortem_interval)] <-
  mean(metadata_filtered$post_mortem_interval, na.rm = TRUE)

metadata_filtered$rin_score[is.na(metadata_filtered$rin_score)] <-
  mean(metadata_filtered$rin_score, na.rm = TRUE)

# Scale covariates.
metadata_filtered$age_at_death <- scale(metadata_filtered$age_at_death)
metadata_filtered$post_mortem_interval <- scale(
  metadata_filtered$post_mortem_interval
)
metadata_filtered$rin_score <- scale(metadata_filtered$rin_score)

# Set reference levels for diagnosis and sex.
metadata_filtered$diagnosis <- factor(
  metadata_filtered$diagnosis,
  levels = c("control", "Alzheimer")
)
metadata_filtered$sex <- factor(
  metadata_filtered$sex,
  levels = c("female", "male")
)


# =============================================================================
# 4. Load gene-count matrix
# =============================================================================

# Gene counts without contaminants.
gene_counts_brain <- read.delim(
  "gene_counts_brain.txt",
  header = TRUE,
  row.names = 1,
  fill = TRUE,
  check.names = FALSE
)


# =============================================================================
# 5. Create and run the DESeq2 model
# =============================================================================

# If a fitted object from a previous run exists, load it; otherwise fit and save.
if (file.exists("dds_dps.rds")) {
  dds_dps <- readRDS("dds_dps.rds")
} else {
  dds_dps <- DESeqDataSetFromMatrix(
    countData = gene_counts_brain,
    colData = metadata_filtered,
    design = ~ age_at_death + rin_score + post_mortem_interval + sex + diagnosis
  )
  
  dds_dps <- DESeq(
    dds_dps,
    test = "Wald",
    fitType = "parametric"
  )
  
  saveRDS(dds_dps, file = "dds_dps.rds")
}


# =============================================================================
# 6. Differential expression: Alzheimer vs control
# =============================================================================

# Check the available result names.
resultsNames(dds_dps)

result_dps_d <- results(
  dds_dps,
  contrast = c("diagnosis", "Alzheimer", "control"),
  pAdjustMethod = "fdr",
  alpha = 0.05
)


# Convert results to a data frame and retain gene IDs.
resultdf_dps_d <- as.data.frame(result_dps_d)
resultdf_dps_d$gene <- rownames(resultdf_dps_d)
resultdf_dps_d <- na.omit(resultdf_dps_d)


# =============================================================================
# 7. Filter differentially expressed genes (DEGs)
# =============================================================================

deg_dps_d <- resultdf_dps_d %>%
  dplyr::filter(
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )

# Optional: save the DEG table.
# write.table(
#   deg_dps_d,
#   "deg_dps_d.txt",
#   quote = FALSE,
#   sep = "\t"
# )


# =============================================================================
# 8. Clean results and count up-/down-regulated genes
# =============================================================================

result_dps_d_df <- as.data.frame(result_dps_d)
result_dps_d_df_clean <- na.omit(result_dps_d_df)

write.table(
  result_dps_d_df_clean,
  "result_dps_d_df_clean.txt",
  quote = FALSE,
  sep = "\t"
)

# Significant genes using adjusted p-value and fold-change thresholds.
result_sig <- result_dps_d_df[
  !is.na(result_dps_d_df$padj) &
    result_dps_d_df$padj < 0.05 &
    abs(result_dps_d_df$log2FoldChange) >= 1,
]

n_up <- sum(result_sig$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(result_sig$log2FoldChange < 0, na.rm = TRUE)

cat("Upregulated genes:", n_up, "\n")
cat("Downregulated genes:", n_down, "\n")
cat("Total:", n_up + n_down, "\n")


# =============================================================================
# 9. MA plot
# =============================================================================

ma_df <- as.data.frame(result_dps_d)
ma_df$Significant <- ifelse(
  ma_df$padj < 0.05 & !is.na(ma_df$padj),
  "FDR < 0.05",
  "NS"
)

p <- ggplot(
  ma_df,
  aes(x = baseMean, y = log2FoldChange, color = Significant)
) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-14, 14)) +
  scale_color_manual(
    values = c("FDR < 0.05" = "#E41A1C", "NS" = "gray")
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "#377eb8"
  ) +
  labs(
    x = "Mean of normalized counts (log10)",
    y = "Log2 fold change",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.text = element_text(size = 11),
    plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm")
  )

ggsave(
  filename = "MA_plot_DESeq2_ggplot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# 10. Volcano plot
# =============================================================================

# Select the top 20 DEGs by adjusted p-value.
top_deg_dps_d <- deg_dps_d %>%
  arrange(padj) %>%
  slice_head(n = 20)

volplot <- EnhancedVolcano(
  resultdf_dps_d,
  lab = resultdf_dps_d$gene,
  selectLab = top_deg_dps_d$gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  max.overlaps = Inf,
  xlab = bquote(~Log[2]~ "Fold Change"),
  ylab = bquote(-Log[10]~ "Adjusted P-value"),
  axisLabSize = 21,
  labSize = 7,
  legendLabSize = 20,
  legendIconSize = 4.0,
  title = "",
  subtitle = NULL,
  caption = NULL,
  legendLabels = c("NS", "Log2FC", "padj", "padj & Log2FC"),
  legendPosition = "right",
  pointSize = 1.2,
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colConnectors = "black",
  arrowheads = FALSE,
  colAlpha = 0.8,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  xlim = c(-48, 48),
  ylim = c(0, 16)
) +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 18)
  )

png(
  "Volcano_plot_diagnosis_with_covariates.png",
  width = 18,
  height = 10,
  units = "in",
  res = 300
)
print(volplot)
dev.off()

cat("SUCCESS: Volcano plot saved.\n")


# =============================================================================
# 11. Heatmap top 50 with ComplexHeatmap
# =============================================================================

#### STEP 1: VST Normalization ####
vsd <- vst(dds_dps, blind = FALSE)
normalized_counts <- assay(vsd)

#### STEP 2: Top 50 DEGs ####
resultdf_dps_d <- na.omit(resultdf_dps_d)
heatmap_top_deg_dps_d <- rownames(deg_dps_d[order(deg_dps_d$padj), ])[1:50]
heatmap_matrix <- normalized_counts[heatmap_top_deg_dps_d, ]

#### STEP 3: Z-score scaling ####
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

#### STEP 4: Sample Annotation ####
sample_annotation <- as.data.frame(colData(dds_dps)[, "sex", drop = FALSE])
colnames(sample_annotation) <- "Sex"
rownames(sample_annotation) <- colnames(heatmap_matrix)

annotation_colors <- list(
  Sex = c(male = "#377eb8", female = "#E41A1C")
)

#### STEP 5: Column Annotation ####
ha_column <- HeatmapAnnotation(
  Sex = sample_annotation$Sex,
  col = annotation_colors,
  show_annotation_name = TRUE,
  annotation_name_side = "left",  # left or right
  annotation_legend_param = list(
    Sex = list(
      title = "Sex",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 10),
      direction = "vertical"
    )
  )
)

#### STEP 6: Heatmap Colors ####
col_fun <- colorRamp2(
  seq(-3, 3, length.out = 50),                  # Z-score limits
  colorRampPalette(c("#2166AC", "white", "#B2182B"))(50)
)

#### STEP 7: Save PNG ####
png(
  filename = "heatmap_top50_dps_d.png",
  width = 3200,
  height = 1800,
  res = 300
)

#### STEP 8: Create Heatmap ####
ht <- Heatmap(
  heatmap_matrix_scaled,
  name = "Z-score",
  col = col_fun,
  
  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  # Row settings
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6),
  
  # Column settings
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 45,
  
  # Annotation
  top_annotation = ha_column,
  
  # Heatmap appearance
  border = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  # Heatmap legend parameters
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9),
    direction = "vertical"
  )
)

#### STEP 9: Draw Heatmap with aligned legends ####
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,                 # Align Sex and Z-score legends
  gap = unit(5, "mm"),                 # Gap between legends
  padding = unit(c(4, 10, 2, 4), "mm") # bottom, right, top, left
)

dev.off()


# =============================================================================
# 12. Session information
# =============================================================================

sessionInfo()

