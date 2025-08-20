required_libraries <- c(
    "data.table",    
    "dplyr",
    "optparse",
    "rlang",
    "ggplot2",
    "purrr",
    "fgsea",
    "RColorBrewer")

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

# Bulk RNA-seq data need to be raw counts.
raw_counts <- fread("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/data/raw_counts_matching_metadata.tsv")
metadata <- fread("/storage/kuijjerarea/ine/breast_met/breast_other/R/data/filtered_metadata_aurora.csv")

