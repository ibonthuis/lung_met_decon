# Functions that are from my single cell - pseudobulk pipeline
meta_per_ct <- function(metadata, celltype_list){
    metadata_ct <- metadata %>%
      dplyr::filter(free_annotation %in% celltype_list) 
   
    return(metadata_ct)
}


filter_count_matrix <- function(count_matrix, metadata_per_ct, celltype_list) {
  # first make the list concatenated again. In that way, you get counts of all cell types of interest in one matrix
  # Then later you need to split the count matrix per celltype again
 # metadata_cts <- purrr::list_rbind(metadata_per_ct)
  idents_to_keep <- metadata_per_ct$V1
  count_matrix <- as.matrix(count_matrix)
  counts <- count_matrix[, colnames(count_matrix) %in% idents_to_keep] 
  return(counts)
}

filter_counts_by_nonzero_geneentries <- function(count_matrix, percent_of_cells) {
    # count_matrix can be of datatype dataframe or matrix
    nr_of_cells <- ncol(count_matrix)*(percent_of_cells/100)
    min_nr_of_cells <- round(nr_of_cells)
    count_matrix <- count_matrix[rowSums(count_matrix != 0) > min_nr_of_cells,]
    cat( ";", paste0(nrow(count_matrix)), " number of genes expressed in ", paste0(percent_of_cells), "% of cells out of a total cell number ", paste0(ncol(count_matrix)), "\n")

    return(count_matrix)
}

filter_counts_by_ct_oi <- function(preprocessed_metadata, oldseurat, celltypes) {
    # preprocessed_metadata$free_annotation <- gsub(" ", "", preprocessed_metadata$free_annotation)
    # metadata_cts_oi <- meta_per_ct(preprocessed_metadata, celltypes)
    # metadata_cts_oi$V1 <- rownames(metadata_cts_oi)
    countmatrix <- oldseurat@raw.data
    countmatrix_oi <- filter_count_matrix(countmatrix, preprocessed_metadata, celltypes)
    return(countmatrix_oi)

}


preprocess_oldseurats_metadata <- function(oldseurat, celltypes_oi){
 # normal_cm <- oldseurat@raw.data
  normal_meta <- oldseurat@meta.data
  normal_meta$free_annotation <- gsub(" ", "", normal_meta$free_annotation)
  normal_meta_cts_oi <- meta_per_ct(normal_meta, celltypes_oi)
  normal_meta_cts_oi$V1 <- rownames(normal_meta_cts_oi)
  return(normal_meta_cts_oi)
  
}

# preprocess_oldseurat_counts <- function(oldseurat, normal_meta, celltypes_oi) {
#   normal_cm <- oldseurat@raw.data
#   normal_counts_cts_oi <- filter_count_matrix(normal_cm, normal_meta, celltypes_oi)


# }