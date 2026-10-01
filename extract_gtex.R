Sys.setenv("VROOM_CONNECTION_SIZE" = 131072000)

library(tidyverse)
library(vroom)

set.seed(42) # For reproducible sampling

# 1. Download & Parse GTEx Sample Attributes
gtex_meta_url <- "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
gtex_meta_file <- "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

if(!file.exists(gtex_meta_file)) {
  message("Downloading GTEx metadata...")
  download.file(gtex_meta_url, gtex_meta_file, mode = "wb")
}

gtex_meta <- read_tsv(gtex_meta_file, show_col_types = FALSE)

# Filter for Whole Blood
blood_samples <- gtex_meta %>%
  filter(SMTSD == "Whole Blood") %>%
  pull(SAMPID)

# Standardize metadata sample IDs just in case they have dots
blood_samples <- gsub("\\.", "-", blood_samples)

# 2. Download GTEx Gene Reads Matrix
gtex_counts_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
gtex_counts_file <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"

if(!file.exists(gtex_counts_file)) {
  message("Downloading GTEx counts matrix (this is ~1.3 GB and may take several minutes)...")
  options(timeout = 3600) # Ensure it doesn't timeout
  download.file(gtex_counts_url, gtex_counts_file, mode = "wb")
}

# 3. Read the entire GTEx counts matrix
message("Reading the full GTEx counts matrix (this might take a moment)...")
gtex_matrix <- read_tsv(
  gtex_counts_file, 
  skip = 2, 
  name_repair = "minimal",
  show_col_types = FALSE
)

# Convert all column names of the loaded GTEx matrix to use dashes instead of dots
colnames(gtex_matrix) <- gsub("\\.", "-", colnames(gtex_matrix))

# Ensure gene column is named 'gene_id'
if ("Name" %in% colnames(gtex_matrix)) {
  gtex_matrix <- gtex_matrix %>% rename(gene_id = Name)
}

# 4. Match and filter the metadata sample IDs against the updated column names
valid_samples <- intersect(blood_samples, colnames(gtex_matrix))

if(length(valid_samples) < 300) {
  stop("Not enough Whole Blood samples matching between metadata and matrix (found ", length(valid_samples), ", need 300)")
}

# Select 300 random samples from the valid intersection
selected_gtex_samples <- sample(valid_samples, 300)
cols_to_keep <- c("gene_id", selected_gtex_samples)

# Subset 300 random samples using select(all_of(cols_to_keep))
message("Subsetting 300 normal samples from GTEx counts matrix...")
gtex_subset <- gtex_matrix %>% select(all_of(cols_to_keep))

# Clean GTEx Gene IDs: Strip version from gene_id
gtex_subset <- gtex_subset %>%
  mutate(gene_id = str_replace(gene_id, "\\..*", "")) %>%
  distinct(gene_id, .keep_all = TRUE)

# 5. Export Final Matrix
out_file <- "gtex_whole_blood_counts.tsv"
message(paste("Exporting to", out_file, "..."))
write_tsv(gtex_subset, out_file)

message(sprintf("Successfully exported final GTEx matrix to %s with %d genes and %d samples.", 
                out_file, nrow(gtex_subset), ncol(gtex_subset)-1))
