# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(sva)
  library(EnhancedVolcano)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Custom theme to remove grid lines as requested
theme_no_grid <- theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

message("Loading data...")
# Load matrix and results
counts_df <- read_tsv("final_300T_300N_matrix.tsv", show_col_types = FALSE)
res_df <- read_tsv("deseq2_myeloma_results.tsv", show_col_types = FALSE)
top_genes <- read_csv("top_significant_myeloma_genes.csv", show_col_types = FALSE)

# Prepare matrix and metadata
count_mat <- as.matrix(counts_df[,-1])
rownames(count_mat) <- counts_df$gene_id

sample_ids <- colnames(count_mat)
condition <- factor(ifelse(grepl("^Myeloma_Tumor_", sample_ids), "Tumor", "Normal"))
batch <- factor(ifelse(grepl("^Myeloma_Tumor_", sample_ids), "MMRF", "GTEx"))

metadata <- data.frame(row.names = sample_ids, condition = condition, batch = batch)

message("Running ComBat_seq to generate corrected counts for visualization...")
# We must remove zero variance genes before ComBat_seq and PCA
variances <- apply(count_mat, 1, var)
count_mat <- count_mat[variances > 0, ]

adjusted_counts <- ComBat_seq(counts = count_mat, batch = metadata$batch)

message("Creating DESeqDataSet and applying Variance Stabilizing Transformation (VST)...")
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts, colData = metadata, design = ~ condition)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# 1. PCA Plot
message("Generating PCA Plot...")
pca <- prcomp(t(vsd_mat))
pca_data <- data.frame(
  PC1 = pca$x[,1], 
  PC2 = pca$x[,2], 
  Condition = metadata$condition
)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(level = 0.95) +
  labs(title = "PCA Plot (VST counts)", x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"), y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_no_grid

ggsave("pca_plot.png", pca_plot, width = 8, height = 6)

# 2. Enhanced Volcano Plot
message("Generating Volcano Plot...")
# We need to map symbols for Volcano
# Filter NA padj
res_volc <- res_df %>% filter(!is.na(padj))
# If symbol is NA, use gene_id
res_volc$symbol[is.na(res_volc$symbol)] <- res_volc$gene_id[is.na(res_volc$symbol)]

volcano_plot <- EnhancedVolcano(res_volc,
  lab = res_volc$symbol,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1.5,
  title = 'Tumor vs Normal',
  subtitle = 'Global DEG Overview',
  pointSize = 1.5,
  labSize = 4,
  selectLab = head(top_genes$symbol, 10),
  drawConnectors = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8)

# 3. Top DEG Hierarchical Heatmap
message("Generating Top 50 DEGs Heatmap...")
top50_ids <- head(top_genes$gene_id, 50)
top50_symbols <- head(top_genes$symbol, 50)
top50_symbols[is.na(top50_symbols)] <- top50_ids[is.na(top50_symbols)]

mat_top50 <- vsd_mat[top50_ids, ]
rownames(mat_top50) <- top50_symbols

annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- rownames(metadata)

png("top50_degs_heatmap.png", width = 1000, height = 800, res = 120)
pheatmap(mat_top50, scale = "row", annotation_col = annotation_col, show_colnames = FALSE, main = "Top 50 Significant DEGs")
dev.off()

# 4. MA Plot
message("Generating MA Plot...")
res_df_ma <- res_df %>% mutate(is_sig = ifelse(padj < 0.05 & abs(log2FoldChange) >= 1.5, "Significant", "Not Significant"))
res_df_ma$is_sig[is.na(res_df_ma$is_sig)] <- "Not Significant"

ma_plot <- ggplot(res_df_ma, aes(x = log10(baseMean+1), y = log2FoldChange, color = is_sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(title = "MA Plot", x = "log10(baseMean + 1)", y = "log2FoldChange") +
  theme_no_grid

ggsave("ma_plot.png", ma_plot, width = 8, height = 6)

# 5. Sample-to-Sample Distance Heatmap
message("Generating Sample Distance Heatmap...")
sampleDists <- dist(t(vsd_mat))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata$condition, rownames(metadata), sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

png("sample_distance_heatmap.png", width = 800, height = 800, res = 120)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors, main = "Sample-to-Sample Distances")
dev.off()

# 6. Top Biomarker Expression Boxplots
message("Generating Top Candidate Genes Boxplots...")
top6_ids <- head(top_genes$gene_id, 6)
top6_symbols <- head(top_genes$symbol, 6)
top6_symbols[is.na(top6_symbols)] <- top6_ids[is.na(top6_symbols)]

top6_mat <- vsd_mat[top6_ids, ]
rownames(top6_mat) <- top6_symbols

top6_df <- as.data.frame(t(top6_mat)) %>% 
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

boxplots <- ggplot(top6_df, aes(x = condition, y = Expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = "Top 6 Candidate Genes Expression", y = "VST Normalized Expression", x = "Condition") +
  theme_no_grid +
  theme(legend.position = "none")

ggsave("top_candidate_genes_boxplots.png", boxplots, width = 10, height = 6)

# 7 & 8. Enrichment Analysis
message("Performing GO and KEGG Enrichment Analysis...")
entrez_ids <- mapIds(org.Hs.eg.db, keys = top_genes$gene_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# GO BP
go_res <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
if (!is.null(go_res) && nrow(go_res) > 0) {
  go_dot <- dotplot(go_res, showCategory = 15, title = "GO BP Enrichment") + theme_no_grid
  ggsave("go_enrichment_dotplot.png", go_dot, width = 10, height = 8)
  
  go_bar <- barplot(go_res, showCategory = 15, title = "GO BP Enrichment") + theme_no_grid
  ggsave("go_enrichment_barplot.png", go_bar, width = 10, height = 8)
} else {
  message("No significant GO terms found.")
}

# KEGG
kegg_res <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.05)
if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  kegg_dot <- dotplot(kegg_res, showCategory = 15, title = "KEGG Pathway Enrichment") + theme_no_grid
  ggsave("kegg_enrichment_dotplot.png", kegg_dot, width = 10, height = 8)
  
  kegg_bar <- barplot(kegg_res, showCategory = 15, title = "KEGG Pathway Enrichment") + theme_no_grid
  ggsave("kegg_enrichment_barplot.png", kegg_bar, width = 10, height = 8)
} else {
  message("No significant KEGG pathways found.")
}

message("All visualizations generated successfully.")