# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(sva)
  library(DESeq2)
  library(org.Hs.eg.db)
})

message("Loading final 300T + 300N count matrix...")
counts_df <- read_tsv("final_300T_300N_matrix.tsv", show_col_types = FALSE)

# Convert to matrix
count_mat <- as.matrix(counts_df[,-1])
rownames(count_mat) <- counts_df$gene_id
message(sprintf("  Matrix: %d genes x %d samples", nrow(count_mat), ncol(count_mat)))

# Build metadata from column names
sample_ids <- colnames(count_mat)
condition  <- ifelse(grepl("^Myeloma_Tumor_", sample_ids), "Myeloma_Tumor", "Normal_Blood")
batch      <- ifelse(grepl("^Myeloma_Tumor_", sample_ids), "MMRF", "GTEx")

metadata <- tibble(
  sample_id = sample_ids,
  condition = condition,
  batch     = batch
)

message(sprintf("  Tumor samples  : %d", sum(metadata$condition == "Myeloma_Tumor")))
message(sprintf("  Normal samples : %d", sum(metadata$condition == "Normal_Blood")))

metadata$condition <- factor(metadata$condition)
metadata$batch     <- factor(metadata$batch)

# Filter out genes with zero variance before PCA
variances <- apply(count_mat, 1, var)
count_mat <- count_mat[variances > 0, ]
message(sprintf("  Genes after removing zero-variance genes: %d", nrow(count_mat)))

# PCA BEFORE Batch Correction (Raw PCA plot)
message("Performing PCA on RAW counts (log2 transformed)...")
log2_counts <- log2(count_mat + 1)

pca_raw <- prcomp(t(log2_counts), scale. = TRUE)
pca_raw_data <- data.frame(
    PC1 = pca_raw$x[,1], 
    PC2 = pca_raw$x[,2], 
    Condition = metadata$condition,
    Batch = metadata$batch,
    Sample = metadata$sample_id
)

pca_raw_plot <- ggplot(pca_raw_data, aes(x = PC1, y = PC2, color = Condition, shape = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = "PCA of Raw Log2 Counts (Before Batch Correction)",
         x = paste0("PC1 (", round(summary(pca_raw)$importance[2,1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_raw)$importance[2,2] * 100, 1), "%)"))

ggsave("pca_plot_raw.png", plot = pca_raw_plot, width = 8, height = 6)
message("Saved Raw PCA plot to pca_plot_raw.png")


# 1. Batch Correction using ComBat_seq
message("Performing batch correction with ComBat_seq...")
adjusted_counts <- ComBat_seq(
  counts = count_mat,
  batch  = metadata$batch
)

# PCA AFTER Batch Correction 
message("Performing PCA on ComBat_seq adjusted counts (log2 transformed)...")
log2_adj_counts <- log2(adjusted_counts + 1)

pca_adj <- prcomp(t(log2_adj_counts), scale. = TRUE)
pca_adj_data <- data.frame(
    PC1 = pca_adj$x[,1], 
    PC2 = pca_adj$x[,2], 
    Condition = metadata$condition,
    Batch = metadata$batch,
    Sample = metadata$sample_id
)

pca_adj_plot <- ggplot(pca_adj_data, aes(x = PC1, y = PC2, color = Condition, shape = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = "PCA of Adjusted Log2 Counts (After ComBat_seq)",
         x = paste0("PC1 (", round(summary(pca_adj)$importance[2,1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_adj)$importance[2,2] * 100, 1), "%)"))

ggsave("pca_plot_adjusted.png", plot = pca_adj_plot, width = 8, height = 6)
message("Saved Adjusted PCA plot to pca_plot_adjusted.png")

# 2. DESeq2 Pipeline
message("Building DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = adjusted_counts,
  colData   = metadata,
  design    = ~ condition
)

# Set Normal_Blood as the reference level
message("Setting 'Normal_Blood' as reference level...")
dds$condition <- relevel(dds$condition, ref = "Normal_Blood")

# Pre-filter: keep genes with at least 10 counts in >= 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message(sprintf("  Genes after DESeq2 pre-filter: %d", nrow(dds)))

message("Running DESeq2 (this may take a while)...")
dds <- DESeq(dds)

# Extract results
res    <- results(dds)
res_df <- as.data.frame(res) %>% rownames_to_column("gene_id")

# 3. Gene Annotation
message("Annotating Ensembl IDs with HUGO symbols...")
res_df$symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys     = res_df$gene_id,
  column   = "SYMBOL",
  keytype  = "ENSEMBL",
  multiVals = "first"
)

# Reorder columns
res_df <- res_df %>%
  dplyr::select(gene_id, symbol, baseMean, log2FoldChange, pvalue, padj, everything())

# 4. Export Full Results
message("Exporting full DESeq2 results...")
write_tsv(res_df, "deseq2_myeloma_results.tsv")

# 5. Export Top DEGs
message("Filtering top significant DEGs (|log2FC| >= 1.5, padj <= 0.05)...")
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(abs(log2FoldChange) >= 1.5, padj <= 0.05) %>%
  arrange(padj)

write_csv(top_genes, "top_significant_myeloma_genes.csv")

message(sprintf("Analysis complete. Found %d significant DEGs.", nrow(top_genes)))
message("Output files:")
message("  - pca_plot_raw.png")
message("  - pca_plot_adjusted.png")
message("  - deseq2_myeloma_results.tsv")
message("  - top_significant_myeloma_genes.csv")