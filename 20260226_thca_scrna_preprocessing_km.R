library(dplyr)
library(Seurat)

# 1. Data Loading and Preprocessing ####
parent_dir <- "thyroid_scrnaseq/"

mt_cutoffs <- c(
  ATC1 = 15, ATC2 = 19, ATC3 = 12,
  PTCEO = 21, PTC1 = 15, PTC2 = 25,
  PTC3 = 25, PTC4 = 20
)

sample_dirs <- list.dirs(parent_dir, recursive = FALSE)

seurat_list <- lapply(sample_dirs, function(dir_path) {
  sample_name <- basename(dir_path)
  counts <- Read10X(data.dir = file.path(dir_path, "filtered_feature_bc_matrix"))
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(
    obj,
    subset = nFeature_RNA > 200 & percent.mt < mt_cutoffs[sample_name]
  )
  return(obj)
})

names(seurat_list) <- names(mt_cutoffs)


# 2. Doublet Detection and Removal ####
library(scDblFinder)
library(DoubletFinder)

doublet_params <- list(
  ATC1 = list(PCs = 1:10, pK = 0.27),
  ATC2 = list(PCs = 1:15, pK = 0.005),
  ATC3 = list(PCs = 1:10, pK = 0.26),
  PTCEO = list(PCs = 1:15, pK = 0.23),
  PTC1 = list(PCs = 1:15, pK = 0.01),
  PTC2 = list(PCs = 1:15, pK = 0.02),
  PTC3 = list(PCs = 1:15, pK = 0.005),
  PTC4 = list(PCs = 1:15, pK = 0.03)
)

run_doublet_detection <- function(obj, PCs, pK) {
  set.seed(123)
  db <- scDblFinder(obj@assays$RNA@counts)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = PCs)
  
  nExp <- round(0.12 * ncol(obj))
  obj <- doubletFinder_v3(
    obj,
    PCs = PCs,
    pN = 0.25,
    pK = pK,
    nExp = nExp,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  
  doublets <- intersect(
    rownames(db@colData)[db@colData$scDblFinder.class == "doublet"],
    rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Doublet"]
  )
  
  obj_qc <- obj[, !colnames(obj) %in% doublets]
  
  return(obj_qc)
}

seurat_list_QC <- mapply(
  function(obj, name) {
    run_doublet_detection(
      obj,
      PCs = doublet_params[[name]]$PCs,
      pK  = doublet_params[[name]]$pK
    )
  },
  seurat_list,
  names(seurat_list),
  SIMPLIFY = FALSE
)


# 3. Integrated Clustering Across Samples ####
THCA_QC <- merge(
  seurat_list_QC[[1]],
  y = seurat_list_QC[-1],
  add.cell.ids = names(seurat_list_QC)
)

THCA_QC <- THCA_QC %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:15)

saveRDS(THCA_QC, "thyroid_scrnaseq/processed.rds")
