library(Seurat)
library(dplyr)
library(Matrix)

# 1. Data Loading and Preprocessing ####
data_dir <- "visium_public/raw"

samples <- list.files(data_dir)

read_visium_sample <- function(sample_name) {
  sample_path <- file.path(data_dir, sample_name)
  
  # Construct full paths for files
  mtx_file <- file.path(sample_path, paste0(sample_name, "_visium_matrix.mtx.gz"))
  features_file <- file.path(sample_path, paste0(sample_name, "_visium_features.tsv.gz"))
  barcodes_file <- file.path(sample_path, paste0(sample_name, "_visium_barcodes.tsv.gz"))
  
  # Read expression matrix
  expr_mat <- ReadMtx(
    mtx = mtx_file,
    features = features_file,
    cells = barcodes_file
  )
  
  return(expr_mat)
}

ST_samples <- lapply(samples, read_visium_sample)
names(ST_samples) <- samples

for (i in names(ST_samples)) {
  ST_samples[[i]] <- CreateSeuratObject(ST_samples[[i]],project = i,
                                        assay = "Spatial")
  ST_samples[[i]]$slice <- i
  ST_samples[[i]]$region <- "thyroid"
  img <- Seurat::Read10X_Image(image.dir = paste(data_dir,"/",i,"/spatial",sep = ""), image.name = 'tissue_hires_image.png')
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = ST_samples[[i]])]
  ST_samples[[i]][['image']] <- img
  DefaultAssay(ST_samples[[i]]) <- "Spatial"
  ST_samples[[i]] <- NormalizeData(ST_samples[[i]])
  ST_samples[[i]] <- FindVariableFeatures(ST_samples[[i]])
  ST_samples[[i]] <- ScaleData(ST_samples[[i]])
  ST_samples[[i]] <- RunPCA(ST_samples[[i]])
  ST_samples[[i]] <- FindNeighbors(ST_samples[[i]])
  ST_samples[[i]] <- FindClusters(ST_samples[[i]], resolution = 0.4)
  ST_samples[[i]] <- RunUMAP(ST_samples[[i]], dims = 1:30)
}

saveRDS(ST_samples,"R_objects/thca_visium_split.rds")

ST_merge <- merge(ST_samples[[1]], y = ST_samples[-1],
                  add.cell.ids = gsub(".*_(.*)-(.*)", "\\1\\2", names(ST_samples)))

saveRDS(ST_merge,"R_objects/thca_visium_merge.rds")


# 2. Deconvolution using RCTD####
library(spacexr)

ST_merge <- readRDS("R_objects/20251107_thca_ST_merge.rds")
ST_merge$sample <- gsub(".*_(.*)-(.*)", "\\1\\2", ST_merge$orig.ident)
ST_merge$sample <- factor(ST_merge$sample, levels = c(paste0('N',1:4),paste0('PTC',1:4),paste0('LPTC',1:4),paste0('ATC',1:4)))
names(ST_merge@images) <- c(paste0('N',1:4),paste0('PTC',1:4),paste0('LPTC',1:4),paste0('ATC',1:4))

ref <- readRDS("R_objects/thca_scrna_ref.rds")
ref <- UpdateSeuratObject(ref)
Idents(ref) <- "merged_anno"

library(SoupX)
library(openxlsx)
ref_celltype_quickmarkers <- quickMarkers(ref@assays$RNA@counts, ref$merged_anno, N=300)
write.xlsx(ref_celltype_quickmarkers,'visium_public/20251110_thca_scrna_celltype_quickmarkers_km.xlsx')
ref_celltype_quickmarkers <- read.xlsx('visium_public/20251110_thca_scrna_celltype_quickmarkers_km.xlsx')
gene_selected <- ref_celltype_quickmarkers$gene %>% unique()

# extract information to pass to the RCTD Reference function
counts <- as.matrix(ref@assays$RNA@counts[gene_selected,])
ref$merged_anno <- gsub('/','.',ref$merged_anno)
cluster <- as.factor(ref$merged_anno)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)


samples <- SplitObject(ST_merge, split.by = "sample")

rctd_results <- lapply(names(samples), function(smp) {
  cat("Running RCTD for sample:", smp, "\n")
  
  obj <- samples[[smp]]
  
  # Extract counts and coordinates
  sp_counts <- as.matrix(obj@assays$Spatial@counts)
  coords <- GetTissueCoordinates(obj, image = smp)
  rownames(coords) <- colnames(sp_counts)
  
  # Create SpatialRNA object
  spatial <- SpatialRNA(coords, sp_counts, colSums(sp_counts))
  
  # Create RCTD object
  myRCTD <- create.RCTD(
    spatialRNA = spatial,
    reference = reference,
    max_cores = 8
  )
  
  # Run RCTD (you can also use doublet_mode = "full" for thorough fitting)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  
  # Save or return results
  saveRDS(myRCTD, file = paste0("visium_public/rctd/RCTD_full_", smp, ".rds"))
  return(myRCTD)
})



