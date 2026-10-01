library(tidyverse)
set.seed(42)
message("Loading GTEx Whole Blood counts...")
gtex <- read_tsv("gtex_whole_blood_counts.tsv", show_col_types = FALSE)
message(sprintf("  GTEx: %d genes x %d samples", nrow(gtex), ncol(gtex) - 1))
message("Loading MMRF labeled expression matrix...")
mmrf_full <- read_tsv("labeled_expression_matrix.tsv", show_col_types = FALSE)
tumor_cols <- colnames(mmrf_full)[grepl("^Myeloma_Tumor_", colnames(mmrf_full))]
message(sprintf("  MMRF: %d tumor samples found", length(tumor_cols)))
if (length(tumor_cols) >= 300) {
  tumor_cols_300 <- sample(tumor_cols, 300)
} else {
  stop(sprintf("Only %d tumor samples; need 300.", length(tumor_cols)))
}
mmrf_300 <- mmrf_full %>% select(gene_id, all_of(tumor_cols_300))
mmrf_300 <- mmrf_300 %>% mutate(gene_id = str_replace(gene_id, "\\..*", "")) %>% distinct(gene_id, .keep_all = TRUE)
gtex      <- gtex      %>% mutate(gene_id = str_replace(gene_id, "\\..*", "")) %>% distinct(gene_id, .keep_all = TRUE)
message("Merging on common genes (inner join)...")
merged <- inner_join(mmrf_300, gtex, by = "gene_id")
n_tumor  <- sum(grepl("^Myeloma_Tumor_", colnames(merged)))
n_normal <- sum(grepl("^GTEX-", colnames(merged)))
message(sprintf("  Merged: %d genes, %d tumor + %d normal", nrow(merged), n_tumor, n_normal))
write_tsv(merged, "final_300T_300N_matrix.tsv")
message("Done! Saved final_300T_300N_matrix.tsv")