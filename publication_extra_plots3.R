suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
  library(RUVSeq)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

message("Loading data...")
counts_df <- read_tsv("final_300T_300N_matrix.tsv", show_col_types = FALSE)
count_mat <- as.matrix(counts_df[, -1])
rownames(count_mat) <- counts_df$gene_id

sample_ids <- colnames(count_mat)
condition <- factor(ifelse(grepl("^Myeloma_Tumor_", sample_ids), "Tumor", "Normal"))
metadata <- data.frame(row.names = sample_ids, condition = condition)

old_res <- read_tsv("limma_myeloma_results.tsv", show_col_types = FALSE)
empirical_controls <- old_res %>% filter(adj.P.Val > 0.5) %>% pull(gene_id)

y_raw <- DGEList(counts = count_mat, group = condition)
keep <- filterByExpr(y_raw)
filtered_counts <- count_mat[keep, ]
empirical_controls <- intersect(empirical_controls, rownames(filtered_counts))

set <- newSeqExpressionSet(as.matrix(filtered_counts), phenoData = data.frame(condition, row.names = sample_ids))
set <- RUVg(set, empirical_controls, k = 1)
metadata$W_1 <- pData(set)$W_1

y <- DGEList(counts = filtered_counts, group = condition)
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + condition + W_1, data = metadata)
colnames(design) <- c("Normal", "Tumor", "W_1")

res_ruv <- read_tsv("ruv_limma_results.tsv", show_col_types = FALSE)
top_genes_ruv <- read_tsv("ruv_significant_degs.tsv", show_col_types = FALSE)

# ============================================================
# PLOT 4: Voom Mean-Variance Trend
# ============================================================
message("Generating Voom Mean-Variance Trend Plot...")
png("voom_meanvariance.png", width = 800, height = 600, res = 120)
v <- voom(y, design, plot = TRUE)
title(main = "Voom Mean-Variance Trend (RUV Design)")
invisible(dev.off())
message("Voom plot done.")

# ============================================================
# PLOT 5: Top DEGs Lollipop Plot (Effect Size)
# ============================================================
message("Generating Top DEGs Lollipop Plot...")

top30 <- top_genes_ruv %>%
  slice_head(n = 30) %>%
  arrange(logFC) %>%
  mutate(
    symbol    = factor(symbol, levels = symbol),
    Direction = ifelse(logFC > 0, "Up-regulated", "Down-regulated")
  )

p_lollipop <- ggplot(top30, aes(x = symbol, y = logFC, color = Direction)) +
  geom_segment(aes(x = symbol, xend = symbol, y = 0, yend = logFC), linewidth = 0.8) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_manual(values = c("Up-regulated" = "#f20f0f", "Down-regulated" = "#0e16ec")) +
  coord_flip() +
  labs(
    title    = "Top 30 DEGs by Log2 Fold Change (RUVSeq)",
    subtitle = "Sorted by effect size (logFC)",
    x        = "Gene Symbol",
    y        = bquote(~Log[2] ~ "Fold Change (Tumor vs Normal)")
  ) +
  theme_pub +
  theme(legend.position = "none")

ggsave("top30_lollipop.png", p_lollipop, width = 9, height = 10)
message("Lollipop plot done.")

# ============================================================
# PLOT 6: GSEA (Gene Set Enrichment Analysis)
# ============================================================
message("Generating GSEA Plot...")

gene_list <- res_ruv %>%
  filter(!is.na(neg_log10_p)) %>%
  mutate(rank_stat = sign(logFC) * neg_log10_p) %>%
  arrange(desc(rank_stat))

entrez_map <- mapIds(org.Hs.eg.db,
  keys      = gene_list$gene_id,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
)
gene_list$entrez <- entrez_map[gene_list$gene_id]
gene_list_clean <- gene_list %>% filter(!is.na(entrez))
ranked_vec <- setNames(gene_list_clean$rank_stat, gene_list_clean$entrez)
ranked_vec <- sort(ranked_vec, decreasing = TRUE)
ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]

set.seed(42)
gsea_res <- gseGO(
  geneList     = ranked_vec,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  message(sprintf("GSEA found %d significant GO:BP terms.", nrow(gsea_res)))

  gsea_dot <- dotplot(gsea_res, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign) +
    theme_pub +
    labs(title = "GSEA GO BP (RUVSeq)", subtitle = "Ranked by signed -log10(p) x direction")
  ggsave("gsea_go_dotplot.png", gsea_dot, width = 14, height = 10)

  gsea_ridge <- ridgeplot(gsea_res, showCategory = 12, fill = "pvalue") +
    labs(title = "GSEA Enrichment Score Distributions", subtitle = "Top 12 GO:BP Pathways") +
    theme_pub
  ggsave("gsea_ridge.png", gsea_ridge, width = 12, height = 10)
} else {
  message("No significant GSEA GO:BP terms found at pvalueCutoff = 0.05.")
}

message("Done! All remaining publication plots generated.")