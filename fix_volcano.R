suppressPackageStartupMessages({
  library(tidyverse)
  library(EnhancedVolcano)
})

theme_no_grid <- theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Fix RUV Volcano
res_ruv <- read_tsv("ruv_limma_results.tsv", show_col_types = FALSE)
res_ruv$adj.P.Val[res_ruv$adj.P.Val == 0] <- .Machine$double.xmin
top_genes_ruv <- res_ruv %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% arrange(adj.P.Val)

volcano_plot_ruv <- EnhancedVolcano(res_ruv,
  lab = res_ruv$symbol,
  x = 'logFC',
  y = 'adj.P.Val',
  pCutoffCol = 'adj.P.Val',
  pCutoff = 0.05,
  FCcutoff = 1.5,
  title = 'Tumor vs Normal (limma-voom + RUVSeq)',
  subtitle = 'Global DEG Overview',
  pointSize = 1.5,
  labSize = 4,
  selectLab = head(top_genes_ruv$symbol, 10),
  drawConnectors = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
ggsave("ruv_volcano_plot.png", volcano_plot_ruv, width = 10, height = 8)

# Fix Old Volcano
old_res <- read_tsv("limma_myeloma_results.tsv", show_col_types = FALSE)
old_res_volc <- old_res %>% filter(!is.na(adj.P.Val))
old_res_volc$symbol[is.na(old_res_volc$symbol)] <- old_res_volc$gene_id[is.na(old_res_volc$symbol)]
old_res_volc$adj.P.Val[old_res_volc$adj.P.Val == 0] <- .Machine$double.xmin
top_old_genes <- old_res_volc %>% filter(adj.P.Val < 0.05, abs(logFC) > 1.5) %>% arrange(adj.P.Val)

volcano_plot_old <- EnhancedVolcano(old_res_volc,
  lab = old_res_volc$symbol,
  x = 'logFC',
  y = 'adj.P.Val',
  pCutoffCol = 'adj.P.Val',
  pCutoff = 0.05,
  FCcutoff = 1.5,
  title = 'Tumor vs Normal (limma-voom) [NO RUV]',
  subtitle = 'Global DEG Overview',
  pointSize = 1.5,
  labSize = 4,
  selectLab = head(top_old_genes$symbol, 10),
  drawConnectors = TRUE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
ggsave("volcano_plot_ylim500.png", volcano_plot_old, width = 10, height = 8)