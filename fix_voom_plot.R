suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
  library(RUVSeq)
})

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    axis.title = element_text(face = "bold")
  )

message("Loading data and rebuilding model...")
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

message("Extracting voom mean-variance trend data...")
# Run voom silently and extract the E (expression) and weights
v <- voom(y, design, plot = FALSE)

# Compute per-gene mean log2-CPM and sqrt(sd) of residuals
log_cpm <- v$E
mean_expr <- rowMeans(log_cpm)
# Use the voom weights to compute approximate sqrt(sigma)
# voom stores precision weights; sqrt(1/w) approximates the standard deviation
sqrt_sigma <- sqrt(1 / rowMeans(v$weights))

voom_df <- data.frame(
  mean_log2cpm = mean_expr,
  sqrt_sd      = sqrt_sigma
)

# Fit a lowess trend line (same as voom internally uses)
lo <- lowess(voom_df$mean_log2cpm, voom_df$sqrt_sd, f = 0.5)
trend_df <- data.frame(x = lo$x, y = lo$y)

message("Generating clean ggplot2 Voom Mean-Variance plot...")
p_voom <- ggplot(voom_df, aes(x = mean_log2cpm, y = sqrt_sd)) +
  geom_point(alpha = 0.15, size = 0.8, color = "#444444") +
  geom_line(data = trend_df, aes(x = x, y = y), color = "#e63946", linewidth = 1.2) +
  labs(
    title    = "Voom Mean-Variance Trend",
    subtitle = "RUV-adjusted design  |  Red line = lowess trend  |  n = 600 samples",
    x        = "Average log\u2082 CPM",
    y        = expression(sqrt(sigma) ~~ "(precision weight proxy)")
  ) +
  theme_pub

ggsave("voom_meanvariance.png", p_voom, width = 8, height = 6, dpi = 150)
message("Done.")