library(tidyverse)
library(purrr)

# Define directories
download_dir <- "downloads"
out_dir <- "."

# Ensure output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# List all STAR gene counts files in the downloads directory
files <- list.files(
  path = download_dir,
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No STAR gene counts files found in 'downloads' directory.")
}

# Function to read and process a single file
process_file <- function(file_path) {
  # Extract sample ID from the folder name (usually the GDC UUID)
  sample_id <- basename(dirname(file_path))
  
  # Read the TSV file, skipping header comments
  df <- read_tsv(file_path, comment = "#", show_col_types = FALSE)
  
  df_processed <- df %>%
    # Filter out STAR summary rows starting with 'N_'
    filter(!str_starts(gene_id, "N_")) %>%
    # Strip Ensembl version numbers
    mutate(gene_id = str_replace(gene_id, "\\..*$", "")) %>%
    # Deduplicate PAR genes
    distinct(gene_id, .keep_all = TRUE) %>%
    # Extract raw unstranded counts and rename to the sample ID
    select(gene_id, !!sym(sample_id) := unstranded)
  
  return(df_processed)
}

# Process all files into a list of dataframes
df_list <- map(files, process_file)

# Merge all samples by gene_id using purrr::reduce(full_join)
combined_matrix <- df_list %>%
  reduce(~ full_join(.x, .y, by = "gene_id"))

# Perform an integer check (%% 1 == 0) on the count columns
int_check <- combined_matrix %>%
  select(-gene_id) %>%
  map_lgl(~ all(is.na(.x) | .x %% 1 == 0))

if (!all(int_check)) {
  warning("Not all values in the combined matrix are integers!")
} else {
  message("Integer check passed: All counts are integers.")
}

# Export the master count matrix
output_file <- file.path(out_dir, "combined_expression_matrix.tsv")
write_tsv(combined_matrix, output_file)
message(paste("Master count matrix saved to", output_file))
