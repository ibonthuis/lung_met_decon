required_libraries <- c(
    "data.table",    
    "dplyr",
    "Seurat"
    #"optparse",
    #"rlang",
    #"ggplot2",
    #"purrr",
    #"fgsea",
   # "RColorBrewer"
   )

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

source("lung_met_decon/pre_processing/pre_process_singlecell_fn.R")

normal_lung <- get(load("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/data/droplet_normal_lung_seurat_ntiss10x.P2.anno.20191002.RC4.Robj"))

part2 <- get(load("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/data/droplet_normal_lung_blood_seurat_ntiss10x.P1.anno.20191002.RC4.Robj"))

part3 <- get(load("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/data/droplet_normal_lung_blood_seurat_ntiss10x.P3.anno.20191002.RC4.Robj"))

# There are 3 seurat objects, old ones.

# I need to take the celltypes I am interested in

# Also the cell types for breast cancer that I already have

# From all these things I need the raw counts

# Create a new seurat object for those counts
# ALTERNATIVE:  not a seurat object because bayesprism takes cell-by-gene raw count matrix as input anyway
# It should be in format dense matrix not a dgcmatrix (just use the as.matrix() function)
cts <- c("Fibrosis-linkedfibroblast", "Myofibroblast", "AirwaySmoothMuscle", "VascularSmoothMuscle", "Fibromyocyte", 'AdventitialFibroblast', "AlveolarFibroblast", "Lipofibroblast", "Pericyte", 'Mesothelial')

#normal_cm <- part2@raw.data
normal_cm <- normal_lung@raw.data
head(normal_cm)
str(normal_cm)
nrow(normal_cm)
ncol(normal_cm)
#filter for sparsity

#normal_meta <- part2@meta.data
normal_meta <- normal_lung@meta.data
normal_meta$free_annotation <- gsub(" ", "", normal_meta$free_annotation)
head(normal_meta)
unique(normal_meta$free_annotation)

unique(normal_meta$free_annotation)

# metadata for the normal_lung object per cell types of interest
normal_meta_cts_oi <- meta_per_ct(normal_meta, cts)
unique(normal_meta_cts_oi$free_annotation)


normal_meta_cts_oi$V1 <- rownames(normal_meta_cts_oi)
normal_counts_cts_oi <- filter_count_matrix(normal_cm, normal_meta_cts_oi, cts)
ncol(normal_counts_cts_oi)

# filter count matrices for the 3 objects
# then merge them, let go of the non-matching genes
# Also probably select a minimum amount of cells of each cell type

table(normal_meta_cts_oi$free_annotation)


normal_cm <- normal_lung@raw.data
normal_meta <- normal_lung@meta.data

normal_meta$free_annotation <- gsub(" ", "", normal_meta$free_annotation)
normal_meta_cts_oi <- meta_per_ct(normal_meta, cts)
normal_meta_cts_oi$V1 <- rownames(normal_meta_cts_oi)
normal_counts_cts_oi <- filter_count_matrix(normal_cm, normal_meta_cts_oi, cts)

norm <- list(normal_meta, normal_counts_cts_oi)
str(norm)
save(
    norm,
    file = "/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/lung_met_decon/pre_processing/normal_cm_and_meta.RData"
)


part3_cm <- part3@raw.data
part3_meta <- part3@meta.data
nrow(part3_meta)
part3_fitlered <- filter_counts_by_ct_oi(part3_meta, part3_cm, cts)
ncol(part3_fitlered) #2571


part2_cm <- part2@raw.data
part2_meta <- part2@meta.data

part2_filtered <- filter_counts_by_ct_oi(part2_meta, part2_cm, cts)
ncol(part2_filtered)


# save counts and metadata, especially 


meta <- preprocess_oldseurats_metadata(part2, cts)
head(meta)
nrow(meta)
nrow(normal_meta)
table(meta$free_annotation)


filtered_counts <- filter_counts_by_ct_oi(meta, part2, cts)

new_meta_and_counts <- list(meta, filtered_counts)

save(
    new_meta_and_counts,
    file = "/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/processed_data/part2_and_meta.RData"
)


meta <- preprocess_oldseurats_metadata(part3, cts)
head(meta)
nrow(meta)
nrow(normal_meta)
table(meta$free_annotation)


filtered_counts <- filter_counts_by_ct_oi(meta, part3, cts)

new_meta_and_counts <- list(meta, filtered_counts)

save(
    new_meta_and_counts,
    file = "/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/processed_data/part3_and_meta.RData"
)

meta <- preprocess_oldseurats_metadata(normal_lung, cts)
head(meta)
nrow(meta)
nrow(normal_meta)
table(meta$free_annotation)


filtered_counts <- filter_counts_by_ct_oi(meta, normal_lung, cts)

new_meta_and_counts <- list(meta, filtered_counts)

save(
    new_meta_and_counts,
    file = "/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/processed_data/part1_and_meta.RData"
)



class(filtered_counts)
#Note: if maximum expression value is <50; CIBERSORTx will assume that data are in log space
# In part3 it looks like there are values higher than 50, so 
max(filtered_counts[, 1])
min(filtered_counts)
quantile(filtered_counts)
head(filtered_counts)
colSums(filtered_counts)
nrow(filtered_counts)

#quantiles <- quantile(filtered_counts$P3_7_TCGTAGATCATCGGAT, probs = c(0.25, 0.5, 0.75))
#quantiles

part1 <- get(load("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/processed_data/part1_and_meta.RData")) # It has 26485 rows
str(part1)

part2 <- get(load("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/processed_data/part2_and_meta.RData"))
part3 <- get(load("/storage/kuijjerarea/ine/projects/BRCA_LUNG_MET/processed_data/part3_and_meta.RData"))

# Does rownames function work for a matrix?

setdiff(rownames(filtered_counts), rownames(part1[[2]]))

rownames(filtered_counts)
 rownames(part1[[2]])
nrow(part1[[2]])

merged_counts <- cbind(part1[[2]], part2[[2]])
merged_counts <- cbind(merged_counts, part3[[2]])

ncol(merged_counts)
ncol(part2[[2]])
ncol(part1[[2]])

tail(rownames(part1[[2]]))
tail(rownames(part2[[2]]))

cts_pt1 <- unique(part1[[1]]$free_annotation)
#cts_pt1 <- unique(meta$free_annotation)
cts_pt2 <- unique(part2[[1]]$free_annotation)
cts_pt3 <- unique(part3[[1]]$free_annotation)

cts_all_parts <- list(cts_pt1, cts_pt2, cts_pt3)
cts_all_parts <- unlist(cts_all_parts)

setdiff(cts_all_parts, cts)
setdiff(cts, cts_all_parts)

