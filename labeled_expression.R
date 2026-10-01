library(tidyverse)

# Define file paths
expr_file <- "combined_expression_matrix.tsv"
meta_file <- "clinical_metadata.tsv"
out_file <- "labeled_expression_matrix.tsv"

# Read the matrices
expr_matrix <- read_tsv(expr_file, show_col_types = FALSE)
metadata <- read_tsv(meta_file, show_col_types = FALSE)

# Get the sample IDs from the expression matrix (excluding 'gene_id')
sample_ids <- colnames(expr_matrix)[colnames(expr_matrix) != "gene_id"]

# Create a mapping from id to sample_type, keeping unique rows
meta_mapping <- metadata %>%
  select(id, sample_type) %>%
  distinct(id, .keep_all = TRUE)

# Create a sequential name for each sample
counters <- list()
new_colnames <- character(length(sample_ids))

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  
  # Find sample type in metadata
  sType <- meta_mapping$sample_type[meta_mapping$id == sid]
  
  if (length(sType) == 0 || is.na(sType)) {
    sType <- "Unknown"
  } else {
    sType <- sType[1]
  }
  
  # Format prefix based on sample condition
  type_prefix <- case_when(
    sType == "Primary_Tumor" ~ "Myeloma_Tumor",
    sType == "Solid_Tissue_Normal" ~ "Myeloma_Normal",
    TRUE ~ paste0("Myeloma_", sType)
  )
  
  # Initialize or increment counter for this prefix
  if (is.null(counters[[type_prefix]])) {
    counters[[type_prefix]] <- 1
  } else {
    counters[[type_prefix]] <- counters[[type_prefix]] + 1
  }
  
  # Assign new sequential name
  new_colnames[i] <- paste0(type_prefix, "_", counters[[type_prefix]])
}

# Apply new column names, preserving 'gene_id'
colnames(expr_matrix) <- c("gene_id", new_colnames)

# Export the updated matrix
write_tsv(expr_matrix, out_file)
message(paste("Labeled expression matrix successfully exported to", out_file))
