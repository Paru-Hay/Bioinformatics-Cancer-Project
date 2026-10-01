suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(RUVSeq)
  library(EnhancedVolcano)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggrepel)
})

theme_no_grid <- theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

message("1. Loading raw matrix and previous limma results...")
counts_df <- read_tsv("final_300T_300N_matrix.tsv", show_col_types = FALSE)
count_mat <- as.matrix(counts_df[, -1])
rownames(count_mat) <- counts_df$gene_id

sample_ids <- colnames(count_mat)
condition <- factor(ifelse(grepl("^Myeloma_Tumor_", sample_ids), "Tumor", "Normal"))
metadata <- data.frame(row.names = sample_ids, condition = condition)

old_res <- read_tsv("limma_myeloma_results.tsv", show_col_types = FALSE)

# ----------------------------------------------------------------------
# FIX OLD VOLCANO PLOT FIRST (with ylim=c(0, 500))
# ----------------------------------------------------------------------
message("Regenerating previous Volcano Plot with ylim=c(0, 500)...")
old_res_volc <- old_res %>% filter(!is.na(adj.P.Val))
old_res_volc$symbol[is.na(old_res_volc$symbol)] <- old_res_volc$gene_id[is.na(old_res_volc$symbol)]
old_res_volc$adj.P.Val[old_res_volc$adj.P.Val == 0] <- .Machine$double.xmin

top_old_genes <- old_res %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>%
  arrange(adj.P.Val)

volcano_plot_old <- EnhancedVolcano(old_res_volc,
  lab = old_res_volc$symbol,
  x = "logFC",
  y = "adj.P.Val",
  pCutoffCol = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1.5,
  title = "Tumor vs Normal (limma-voom) [NO RUV]",
  subtitle = "Global DEG Overview",
  pointSize = 1.5,
  labSize = 4,
  selectLab = head(top_old_genes$symbol, 10),
  drawConnectors = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
ggsave("volcano_plot_ylim500.png", volcano_plot_old, width = 10, height = 8)

# ----------------------------------------------------------------------
# NEW RUVSEQ PIPELINE
# ----------------------------------------------------------------------
message("2. Identifying empirical control genes...")
empirical_controls <- old_res %>%
  filter(adj.P.Val > 0.5) %>%
  pull(gene_id)
message(sprintf("Found %d empirical control genes.", length(empirical_controls)))

message("3. Running RUVg (k=1) on raw counts...")
y <- DGEList(counts = count_mat, group = condition)
keep <- filterByExpr(y)
filtered_counts <- count_mat[keep, ]
empirical_controls <- intersect(empirical_controls, rownames(filtered_counts))

set <- newSeqExpressionSet(as.matrix(filtered_counts), phenoData = data.frame(condition, row.names = sample_ids))
set <- RUVg(set, empirical_controls, k = 1)
metadata$W_1 <- pData(set)$W_1

message("4. Running limma-voom DE Analysis with RUV factors...")
y <- DGEList(counts = filtered_counts, group = condition)
y <- calcNormFactors(y)

design <- model.matrix(~ 0 + condition + W_1, data = metadata)
colnames(design) <- c("Normal", "Tumor", "W_1")

v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)
contr <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit.cont <- contrasts.fit(fit, contr)
fit.cont <- eBayes(fit.cont)

df_total <- fit.cont$df.prior + fit.cont$df.residual
neg_log10_p_all <- -(pt(abs(fit.cont$t[, 1]), df = df_total, lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)

res_ruv <- topTable(fit.cont, coef = 1, number = Inf, adjust.method = "BH") %>%
  rownames_to_column("gene_id")

res_ruv$symbol <- mapIds(org.Hs.eg.db, keys = res_ruv$gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_ruv$symbol[is.na(res_ruv$symbol)] <- res_ruv$gene_id[is.na(res_ruv$symbol)]

res_ruv$neg_log10_p <- neg_log10_p_all[res_ruv$gene_id]

write_tsv(res_ruv, "ruv_limma_results.tsv")

top_genes_ruv <- res_ruv %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>%
  arrange(adj.P.Val)
write_tsv(top_genes_ruv, "ruv_significant_degs.tsv")
message(sprintf("RUVSeq Pipeline found %d significant DEGs.", nrow(top_genes_ruv)))

# ----------------------------------------------------------------------
# VISUALIZATIONS FOR RUVSEQ
# ----------------------------------------------------------------------
message("5. Generating DESeq2 VST and removing W_1 effect for visualization...")
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = ~ condition + W_1)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
vsd_mat_adj <- removeBatchEffect(vsd_mat, covariates = metadata$W_1, design = model.matrix(~ metadata$condition))

message("Generating RUV PCA Plot...")
pca <- prcomp(t(vsd_mat_adj))
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = metadata$condition)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(level = 0.95) +
  labs(title = "PCA Plot (RUV adjusted VST counts)", x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"), y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")) +
  theme_no_grid
ggsave("ruv_pca_plot.png", pca_plot, width = 8, height = 6)

message("Generating RUV Volcano Plot...")
res_ruv <- res_ruv %>%
  mutate(sig_category = case_when(
    adj.P.Val < 0.05 & logFC > 1.5 ~ "Up-regulated",
    adj.P.Val < 0.05 & logFC < -1.5 ~ "Down-regulated",
    adj.P.Val < 0.05 & abs(logFC) <= 1.5 ~ "Significant p-value",
    TRUE ~ "NS"
  ))

volcano_plot_ruv <- ggplot(res_ruv, aes(x = logFC, y = neg_log10_p, color = sig_category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Up-regulated" = "#f20f0f",
    "Down-regulated" = "#0e16ec",
    "Significant p-value" = "#A9A9A9",
    "NS" = "#D3D3D3"
  )) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
  geom_text_repel(data = head(top_genes_ruv, 10), aes(label = symbol), color = "black", size = 4, max.overlaps = Inf) +
  labs(
    title = "Tumor vs Normal (limma-voom + RUVSeq)",
    subtitle = "Global DEG Overview",
    x = bquote(~ Log[2] ~ "fold change"),
    y = bquote(~ -Log[10] ~ "p-value (calculated from t-statistic)")
  ) +
  theme_no_grid +
  theme(legend.position = "none") +
  ylim(0, max(res_ruv$neg_log10_p, na.rm = TRUE) + 20)

ggsave("ruv_volcano_plot.png", volcano_plot_ruv, width = 10, height = 8)

message("Generating RUV Top 50 DEGs Heatmap...")
top50_ids <- head(top_genes_ruv$gene_id, 50)
top50_symbols <- head(top_genes_ruv$symbol, 50)
mat_top50 <- vsd_mat_adj[top50_ids, ]
rownames(mat_top50) <- top50_symbols
annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- rownames(metadata)

png("ruv_top50_degs_heatmap.png", width = 1000, height = 800, res = 120)
pheatmap(mat_top50, scale = "row", annotation_col = annotation_col, show_colnames = FALSE, main = "Top 50 DEGs (RUV Adjusted)")
invisible(dev.off())

message("Generating RUV MA Plot...")
res_ruv_ma <- res_ruv %>% mutate(is_sig = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 1.5, "Significant", "Not Significant"))
ma_plot <- ggplot(res_ruv_ma, aes(x = AveExpr, y = logFC, color = is_sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  labs(title = "MA Plot (RUVSeq)", x = "Average Expression (log2-CPM)", y = "log2FoldChange") +
  theme_no_grid
ggsave("ruv_ma_plot.png", ma_plot, width = 8, height = 6)

message("Generating RUV Sample Distance Heatmap...")
sampleDists <- dist(t(vsd_mat_adj))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata$condition, rownames(metadata), sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

png("ruv_sample_distance_heatmap.png", width = 800, height = 800, res = 120)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors, main = "Sample-to-Sample Distances (RUV Adjusted)")
invisible(dev.off())

message("Generating RUV Top Candidate Genes Boxplots...")
top6_ids <- head(top_genes_ruv$gene_id, 6)
top6_symbols <- head(top_genes_ruv$symbol, 6)
top6_mat <- vsd_mat_adj[top6_ids, ]
rownames(top6_mat) <- top6_symbols

top6_df <- as.data.frame(t(top6_mat)) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

boxplots <- ggplot(top6_df, aes(x = condition, y = Expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~Gene, scales = "free_y") +
  labs(title = "Top 6 Candidate Genes Expression (RUV Adjusted)", y = "VST Normalized Expression", x = "Condition") +
  theme_no_grid +
  theme(legend.position = "none")
ggsave("ruv_top_candidate_genes_boxplots.png", boxplots, width = 10, height = 6)

message("Performing Enrichment Analysis on RUV DEGs...")
entrez_ids <- mapIds(org.Hs.eg.db, keys = top_genes_ruv$gene_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# GO
go_res <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
if (!is.null(go_res) && nrow(go_res) > 0) {
  go_dot <- dotplot(go_res, showCategory = 15, title = "GO BP Enrichment (RUVSeq)") + theme_no_grid
  ggsave("ruv_go_enrichment_dotplot.png", go_dot, width = 10, height = 8)

  go_bar <- barplot(go_res, showCategory = 15, title = "GO BP Enrichment (RUVSeq)") + theme_no_grid
  ggsave("ruv_go_enrichment_barplot.png", go_bar, width = 10, height = 8)
}

# KEGG
kegg_res <- enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  kegg_dot <- dotplot(kegg_res, showCategory = 15, title = "KEGG Pathway Enrichment (RUVSeq)") + theme_no_grid
  ggsave("ruv_kegg_enrichment_dotplot.png", kegg_dot, width = 10, height = 8)

  kegg_bar <- barplot(kegg_res, showCategory = 15, title = "KEGG Pathway Enrichment (RUVSeq)") + theme_no_grid
  ggsave("ruv_kegg_enrichment_barplot.png", kegg_bar, width = 10, height = 8)
}

message("Done! All RUVSeq tasks completed.")
