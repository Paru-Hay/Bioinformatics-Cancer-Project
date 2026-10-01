library(tidyverse)

# 1. Load Data
message("Loading count matrices...")
# We use the processed MMRF count matrix
mmrf_file <- "labeled_expression_matrix.tsv" 
if (!file.exists(mmrf_file)) {
  # Fallback to combined if labeled doesn't exist
  mmrf_file <- "combined_expression_matrix.tsv"
}
mmrf_counts <- read_tsv(mmrf_file, show_col_types = FALSE)

gtex_file <- "gtex_whole_blood_counts.tsv"
gtex_counts <- read_tsv(gtex_file, show_col_types = FALSE)

# 2. Feature Alignment
message("Aligning features...")
# Strip decimal extensions from Ensembl IDs and keep unique
mmrf_counts <- mmrf_counts %>%
  mutate(gene_id = str_replace(gene_id, "\\..*$", "")) %>%
  distinct(gene_id, .keep_all = TRUE)

gtex_counts <- gtex_counts %>%
  mutate(gene_id = str_replace(gene_id, "\\..*$", "")) %>%
  distinct(gene_id, .keep_all = TRUE)

# Keep only intersecting genes
common_genes <- intersect(mmrf_counts$gene_id, gtex_counts$gene_id)
message(sprintf("Found %d intersecting genes between datasets.", length(common_genes)))

mmrf_counts <- mmrf_counts %>% filter(gene_id %in% common_genes)
gtex_counts <- gtex_counts %>% filter(gene_id %in% common_genes)

# 3. Matrix Assembly
message("Merging matrices...")
combined_matrix <- inner_join(mmrf_counts, gtex_counts, by = "gene_id")

# 4. Sample Metadata Generation
message("Generating sample metadata...")
mmrf_samples <- colnames(mmrf_counts)[colnames(mmrf_counts) != "gene_id"]
gtex_samples <- colnames(gtex_counts)[colnames(gtex_counts) != "gene_id"]

# Categorize condition (assuming MMRF are tumors and GTEx are normal blood)
metadata <- tibble(
  sample_id = c(mmrf_samples, gtex_samples),
  batch = c(rep("MMRF", length(mmrf_samples)), rep("GTEx", length(gtex_samples)))
) %>%
  mutate(condition = ifelse(batch == "GTEx", "Normal_Blood", "Myeloma_Tumor"))

# 5. Quality Control & Integer Validation
message("Performing quality control...")
# Extract just the count data for QC
count_data <- combined_matrix %>% select(-gene_id)

# Strict integer check: all values must be integers
is_integer <- all(count_data %% 1 == 0, na.rm = TRUE)
if (!is_integer) {
  warning("QC Warning: Not all values in the combined matrix are integers!")
} else {
  message("QC Success: Integer check passed. All counts are integers.")
}

# Filter out low-expression features: >= 10 counts in >= 10 samples
meets_threshold <- count_data >= 10
samples_meeting_threshold <- rowSums(meets_threshold, na.rm = TRUE)
keep_genes <- samples_meeting_threshold >= 10

combined_matrix_filtered <- combined_matrix[keep_genes, ]
message(sprintf("Low-expression filter: Retained %d out of %d genes.", nrow(combined_matrix_filtered), nrow(combined_matrix)))

# 6. Export
message("Exporting final matrix and metadata...")
write_tsv(combined_matrix_filtered, "combined_myeloma_expression_matrix.tsv")
write_tsv(metadata, "myeloma_metadata.tsv")

message("Process completed successfully.")
