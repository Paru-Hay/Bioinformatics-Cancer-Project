suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggridges)
})

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9, lineheight = 0.85),
    strip.text = element_text(face = "bold", size = 11)
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

# --- Wrap long pathway description labels ---
gsea_res2 <- gsea_res
gsea_res2@result$Description <- str_wrap(gsea_res2@result$Description, width = 40)

# Dotplot: show top 12 per direction, increase height, reduce font
gsea_dot <- dotplot(gsea_res2, showCategory = 12, split = ".sign") +
  facet_grid(. ~ .sign) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(
    title    = "GSEA GO:BP Enrichment (RUVSeq)",
    subtitle = "Ranked by signed -log10(p) x fold-change direction"
  ) +
  theme_pub +
  theme(
    axis.text.y   = element_text(size = 8, lineheight = 0.9),
    plot.margin   = margin(10, 20, 10, 10)
  )

ggsave("gsea_go_dotplot.png", gsea_dot,
       width = 16, height = 12, dpi = 150, limitsize = FALSE)
message("GSEA dotplot done.")