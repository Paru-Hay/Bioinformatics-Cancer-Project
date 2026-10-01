library(tidyverse)

message("Loading final matrix to rename columns...")
df <- read_tsv("final_300T_300N_matrix.tsv", show_col_types = FALSE)

cols <- colnames(df)
tumor_idx <- grep("^Myeloma_Tumor_", cols)
normal_idx <- grep("^GTEX-", cols)

if(length(tumor_idx) > 0 && length(normal_idx) > 0) {
  cols[tumor_idx] <- paste0("Myeloma_Tumor_", seq_along(tumor_idx))
  cols[normal_idx] <- paste0("Normal_Sample_", seq_along(normal_idx))
  colnames(df) <- cols
  message("Renaming successful. Writing back to final_300T_300N_matrix.tsv...")
  write_tsv(df, "final_300T_300N_matrix.tsv")
  message("Done.")
} else {
  message("Columns might already be renamed or not found.")
}