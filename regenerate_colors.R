suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

theme_no_grid <- theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

res_ruv <- read_tsv("ruv_limma_results.tsv", show_col_types = FALSE)
top_genes_ruv <- read_tsv("ruv_significant_degs.tsv", show_col_types = FALSE)

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
    "Up-regulated" = "#FF1493",
    "Down-regulated" = "#00CED1",
    "Significant p-value" = "#A9A9A9",
    "NS" = "#D3D3D3"
  )) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
  geom_text_repel(data = head(top_genes_ruv, 10), aes(label = symbol), color = "black", size = 4, max.overlaps = Inf) +
  labs(
    title = 'Tumor vs Normal (limma-voom + RUVSeq)',
    subtitle = 'Global DEG Overview',
    x = bquote(~Log[2]~ 'fold change'),
    y = bquote(~-Log[10]~ 'p-value (calculated from t-statistic)')
  ) +
  theme_no_grid +
  theme(legend.position = "none") +
  ylim(0, max(res_ruv$neg_log10_p, na.rm=TRUE) + 20)

ggsave("ruv_volcano_plot.png", volcano_plot_ruv, width = 10, height = 8)