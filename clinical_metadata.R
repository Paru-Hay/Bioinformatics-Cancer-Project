library(dplyr)
library(readr)
library(stringr)
library(TCGAutils)

# Read the manifest
manifest <- read_tsv("gdc_manifest_myeloma.txt", show_col_types = FALSE)

# Convert UUIDs to barcodes using TCGAutils
# We extract the file_id column which corresponds to the id column in the manifest
mapping <- UUIDtoBarcode(manifest$id, from_type = "file_id")

# Create a mapping dictionary or join it
# The output of UUIDtoBarcode is a data.frame with columns file_id and aliquot_submitter_id (which is the barcode)
# We will join it with the manifest on file_id == id
mapped_manifest <- manifest %>%
  left_join(mapping, by = c("id" = "file_id")) %>%
  rename(barcode = associated_entities.entity_submitter_id)

# Extract characters 14-15 of barcode to identify sample types
mapped_manifest <- mapped_manifest %>%
  mutate(
    sample_code = str_sub(barcode, 14, 15),
    sample_type = case_when(
      sample_code == "M_" ~ "Primary_Tumor",
      sample_code == "B_" ~ "Solid_Tissue_Normal",
      TRUE ~ "Other"
    )
  )

# Save the output as TSV
write_tsv(mapped_manifest, "clinical_metadata.tsv")
