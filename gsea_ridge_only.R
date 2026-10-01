suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggridges)
})

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

message("Building ranked gene list...")
res_ruv <- read_tsv("ruv_limma_results.tsv", show_col_types = FALSE)

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

message("Running GSEA...")
set.seed(42)
gsea_res <- gseGO(
  geneList     = ranked_vec,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  eps          = 0
)

message(sprintf("GSEA found %d significant GO:BP terms.", nrow(gsea_res)))

gsea_ridge <- ridgeplot(gsea_res, showCategory = 12, fill = "pvalue") +
  labs(title = "GSEA Enrichment Score Distributions", subtitle = "Top 12 GO:BP Pathways") +
  theme_pub
ggsave("gsea_ridge.png", gsea_ridge, width = 12, height = 10)

message("Ridge plot done.")