#Fig. 3A - TF activity ####
library(dplyr)
library(Seurat)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)

epi <- readRDS('R_objects/THCA_scrna_epithelial.rds')
exprMat <- as.matrix(epi@assays$RNA@counts)
saveRDS(exprMat, "int/exprMat.rds")

#initialize SCENIC settings
org <- "hgnc"
dbDir <- "/home/jc2545/palmer_scratch/km/thca/SCENIC/cisTarget_databases/"
myDatasetTitle <- "thca_epi_sub3"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#scenic_corr_genie3
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)
runCorrelation(exprMat_filtered, scenicOptions)
dim(exprMat_filtered)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

#build and score the GRN
readRDS("./exprMat.rds") -> exprMat
readRDS("./int/scenicOptions.Rds") -> scenicOptions
exprMat_log <- log2(exprMat+1)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# downstream analysis
regulonAUC <- readRDS('scenic_epi_sub3_macbook/int/3.4_regulonAUC.Rds')

cancertypes <- epi@meta.data %>% select(annotation_epi)
cancertypes$annotation_epi <- factor(cancertypes$annotation_epi, levels = c("Follicular cell", "PTC type1", "PTC type2", "ATC"))
selectedResolution <- "annotation_epi"
cellsPerTypes <- split(rownames(cancertypes), cancertypes[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_bycancertypes <- sapply(cellsPerTypes,
                                        function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_bycancertypes_Scaled <- t(scale(t(regulonActivity_bycancertypes), center = T, scale=T))

hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_bycancertypes_Scaled, name="Regulon activity", cluster_columns = F,
                                   clustering_method_rows = 'average',
                                   row_dend_width = unit(4, "cm"),
                                   row_names_gp=grid::gpar(fontsize=8),
                                   column_names_gp=grid::gpar(fontsize=16),
                                   column_names_rot = 45))

# Filtering regulons highly activated in ATC or PTC
#Follicular cell: 240 - 258
#PTC type1: 216 - 221
#PTC type2: 50 - 202
#ATC: 12 - 28
interest_regulons <- row.names(regulonActivity_bycancertypes_Scaled)[row_order(hm)[c(240:258,216:221,50:202,12:28)]]

hm2 <- draw(ComplexHeatmap::Heatmap(regulonActivity_bycancertypes_Scaled[interest_regulons,], name="Regulon activity",
                                    show_column_dend = F, cluster_columns = F, clustering_method_rows = 'average',
                                    row_names_gp=grid::gpar(fontsize=14, fontface = "bold"),
                                    row_dend_width = unit(65, "mm"),
                                    column_names_rot = 45,
                                    column_names_gp = grid::gpar(fontsize=20, fontface = "bold"),
                                    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, 'cm'),
                                                                title_gp = grid::gpar(fontsize=20, fontface = "bold"),
                                                                labels_gp = grid::gpar(fontsize=15)),
                                    width = unit(20, "cm"), height = unit(80, "cm")
), heatmap_legend_side = "top")
clustered_interest_regulons <- row.names(regulonActivity_bycancertypes_Scaled[interest_regulons,])[row_order(hm2)]

# Filtering regulons highly activated in ATC or PTC using wilcoxon test
regulonAUC_scaled <- t(scale(t(getAUC(regulonAUC)), center = T, scale=T))
cellorder <- c(cellsPerTypes$`Follicular cell`, cellsPerTypes$`PTC type1`, 
               cellsPerTypes$`PTC type2`, cellsPerTypes$ATC)
cluster_anno <- data.frame(cluster = rep(c('Follicular', 'PTC type1', 'PTC type2', 'ATC'), c(847,1230,1957,1905)))
regulonAUC_scaled_order <- regulonAUC_scaled[clustered_interest_regulons,cellorder]
row.names(cluster_anno) <- colnames(regulonAUC_scaled_order)

cluster_sample_anno_colors = list(cluster = c(Follicular = "#EA8331", `PTC type1` = '#7CAE00', 
                                              `PTC type2` = "#35A2FF", ATC = "#C00000"),
                                  sample = c(PTC2 = '#F8766D', PTC3 = '#C09B00', PTC4 = '#00BF7D', 
                                             ATC1 = '#00BAE0', ATC2 = '#9590FF', ATC3 = '#E76BF3')) #scales::show_col(scales::hue_pal()(10))

#regulons significantly activated in Follicular cells
clustered_interest_regulons_nameonly <- as.character(do.call(rbind.data.frame, strsplit(clustered_interest_regulons, " "))[[1]])
wilcox_follicular <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("row_number","Regulon", "p.value"))
for (i in 1:19){
  res <- wilcox.test(regulonAUC_scaled_order[i,cellsPerTypes$`Follicular cell`], regulonAUC_scaled_order[i,setdiff(cellorder, cellsPerTypes$`Follicular cell`)])
  if(res$p.value < 0.05) {
    print(paste0(i, ': ', clustered_interest_regulons_nameonly[i]))
    res_v2 <- data.frame(row_number = i, Regulon = clustered_interest_regulons_nameonly[i], p.value = res$p.value)
    wilcox_follicular <- rbind(wilcox_follicular, res_v2)
  }
}

#regulons significantly activated in PTC type1
clustered_interest_regulons_nameonly <- as.character(do.call(rbind.data.frame, strsplit(clustered_interest_regulons, " "))[[1]])
wilcox_PTCtype1 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("row_number","Regulon", "p.value"))
for (i in 173:178){
  res <- wilcox.test(regulonAUC_scaled_order[i,cellsPerTypes$`PTC type1`], regulonAUC_scaled_order[i,setdiff(cellorder, cellsPerTypes$`PTC type1`)])
  if(res$p.value < 0.05) {
    print(paste0(i, ': ', clustered_interest_regulons_nameonly[i]))
    res_v2 <- data.frame(row_number = i, Regulon = clustered_interest_regulons_nameonly[i], p.value = res$p.value)
    wilcox_PTCtype1 <- rbind(wilcox_PTCtype1, res_v2)
  }
}

#regulons significantly activated in PTC type2
clustered_interest_regulons_nameonly <- as.character(do.call(rbind.data.frame, strsplit(clustered_interest_regulons, " "))[[1]])
wilcox_PTCtype2 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("row_number","Regulon", "p.value"))
for (i in 20:172){
  res <- wilcox.test(regulonAUC_scaled_order[i,cellsPerTypes$`PTC type2`], regulonAUC_scaled_order[i,setdiff(cellorder, cellsPerTypes$`PTC type2`)])
  if(res$p.value < 0.05) {
    print(paste0(i, ': ', clustered_interest_regulons_nameonly[i]))
    res_v2 <- data.frame(row_number = i, Regulon = clustered_interest_regulons_nameonly[i], p.value = res$p.value)
    wilcox_PTCtype2 <- rbind(wilcox_PTCtype2, res_v2)
  }
}

#regulons significantly activated in ATC
clustered_interest_regulons_nameonly <- as.character(do.call(rbind.data.frame, strsplit(clustered_interest_regulons, " "))[[1]])
wilcox_ATChigh <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("row_number","Regulon", "p.value"))
for (i in 179:195){
  res <- wilcox.test(regulonAUC_scaled_order[i,cellsPerTypes$ATC], regulonAUC_scaled_order[i,setdiff(cellorder, cellsPerTypes$ATC)])
  if(res$p.value < 0.05) {
    print(paste0(i, ': ', clustered_interest_regulons_nameonly[i]))
    res_v2 <- data.frame(row_number = i, Regulon = clustered_interest_regulons_nameonly[i], p.value = res$p.value)
    wilcox_ATChigh <- rbind(wilcox_ATChigh, res_v2)
  }
}


wilcox_all <- do.call("rbind",list(wilcox_follicular,wilcox_PTCtype1, wilcox_PTCtype2, wilcox_ATChigh))

#significantly activated regulons (p.value < 1e-220)
wilcox_all_sig <- wilcox_all %>% filter(p.value < 1e-220)

interest_regulons_significant_row_number <-  wilcox_all_sig$row_number

interest_regulons_significant <- row.names(regulonAUC_scaled_order)[interest_regulons_significant_row_number]

hm3 <- draw(ComplexHeatmap::Heatmap(regulonActivity_bycancertypes_Scaled[interest_regulons_significant,], name="Regulon activity",
                                    show_column_dend = F, cluster_columns = F, clustering_method_rows = 'mcquitty',
                                    row_names_gp=grid::gpar(fontsize=25, fontface = "bold"),
                                    row_dend_width = unit(65, "mm"),
                                    column_names_rot = 45,
                                    column_names_gp = grid::gpar(fontsize=23, fontface = "bold"),
                                    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, 'cm'),
                                                                title_gp = grid::gpar(fontsize=18, fontface = "bold"),
                                                                labels_gp = grid::gpar(fontsize=15)),
                                    width = unit(20, "cm"), height = unit(80, "cm")
), heatmap_legend_side = "top")

interest_regulons_significant_clustered <- row.names(regulonActivity_bycancertypes_Scaled[interest_regulons_significant,])[row_order(hm3)]
write.csv(interest_regulons_significant_clustered, 'scenic_epi_sub3_macbook/exploring_data_including_follicular/thca_by_cancertype_clustering_average/20240417_thca_by_cancertype_interesting_clustering_wilcox_mcquitty_regulons_km.csv')

# TF activity by cell unit
interest_regulons_significant_clustered <- read.csv('scenic_epi_sub3_macbook/exploring_data_including_follicular/thca_by_cancertype_clustering_average/20240417_thca_by_cancertype_interesting_clustering_wilcox_mcquitty_regulons_km.csv',
                                                    row.names = 1)
regulonAUC_scaled_order_wilcox <- regulonAUC_scaled_order[interest_regulons_significant_clustered$x,]

reordered_cells <- c(sample(cellsPerTypes$`Follicular cell`),
                     sample(cellsPerTypes$`PTC type1`),
                     sample(cellsPerTypes$`PTC type2`),
                     sample(cellsPerTypes$ATC))

pheatmap::pheatmap(regulonAUC_scaled_order_wilcox[,reordered_cells], border_color = NA, color = my_color, breaks = seq(-1.5,1.5, length.out=100), 
                   scale = 'none', cluster_cols = F, clustering_method = 'ward.D2', show_colnames = F,
                   annotation_col = cluster_anno, annotation_colors = cluster_sample_anno_colors, 
                   fontsize = 10, fontsize_row = 12)


#Fig. 3B - pathways enriched in CEBPB regulon ####
library(enrichR)
library(ggplot2)
library(viridis)

regulons <- readRDS('int/2.6_regulons_asGeneSet.Rds')
cebpb <- regulons$CEBPB

dbs <- c("MSigDB_Hallmark_2020")
df_enriched <- as.data.frame(enrichr(cebpb, dbs))
colnames(df_enriched) <- gsub("MSigDB_Hallmark_2020.", "", colnames(df_enriched))
df_enriched$Genes <- gsub(";", ", ", df_enriched$Genes)
df_enriched_flt <- df_enriched %>% select(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>% filter(Adjusted.P.value < 0.05) %>% arrange(desc(Combined.Score))

count <- as.numeric(do.call(rbind.data.frame, strsplit(df_enriched_flt$Overlap, "/"))[[1]])
total <- as.numeric(do.call(rbind.data.frame, strsplit(df_enriched_flt$Overlap, "/"))[[2]])
df_enriched_flt$Gene.Ratio <- count/total
df_enriched_flt$Count <- count
colnames(df_enriched_flt)[3] <- 'Adjusted P.value'

ggplot(df_enriched_flt, aes(x = Gene.Ratio, y = reorder(Term, Gene.Ratio), color = `Adjusted P.value`, size = Count)) +
  geom_point() +
  scale_size(range = c(1,6), breaks = seq(0,30,5)) +
  scale_color_viridis(direction = -1) +
  theme_linedraw() +
  labs(y = "", x = "Gene Ratio") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "top", 
        legend.text = element_text(size = 16.5),
        legend.title = element_text(size = 16.5), 
        axis.title.x = element_text(size = 16.5),
        axis.text.x = element_text(size = 16.5, color = "black"),
        axis.ticks = element_line(color = "black", size=2))


#Fig. 3C - correlation between IL-6/JAK/STAT signature, EMT, and ATC signature ####
library(ggpubr)
#IL-6/JAK/STAT - EMT
emt_cebpb_regulon <- df_enriched %>% filter(Term == 'Epithelial Mesenchymal Transition') %>% pull(Genes) %>% strsplit(., ", ")
jakstat_cebpb_regulon <- df_enriched %>% filter(Term == 'IL-6/JAK/STAT3 Signaling') %>% pull(Genes) %>% strsplit(., ", ")
epi <- AddModuleScore(epi, features = list(emt_cebpb_regulon[[1]]), name="emt_cebpb_regulon")
epi <- AddModuleScore(epi, features = list(jakstat_cebpb_regulon[[1]]), name="jakstat_cebpb_regulon")

jakstat_emt <- data.frame(
  jakstat = epi$jakstat_cebpb_regulon1,
  emt = epi$emt_cebpb_regulon1,
  group = epi$annotation_epi)
jakstat_emt$group <- factor(jakstat_emt$group, levels = c('Follicular cell', 'PTC type1', 'PTC type2', 'ATC'))

result <- cor.test(jakstat_emt$jakstat, jakstat_emt$emt)

ggplot(jakstat_emt, aes(x = jakstat, y = emt)) +
  geom_point(aes(shape = group, color = group)) +
  scale_color_manual(values = c('#EA8331', '#7CAE00', '#35A2FF', '#C00000')) +
  geom_smooth(method = 'lm', color ='red') +
  annotate("text", x = 0.9, y = 1.35, col = "black", size = 10,
           label = paste("r = ", signif(cor(jakstat_emt$jakstat, jakstat_emt$emt),3)))+
  theme_classic() +
  labs(x = 'JAK STAT Signature', y = 'EMT Signature') +
  theme(axis.title.x = element_text(size = 24, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 24, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 24, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 24, family = 'Arial', colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 24, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank()) &
  guides(shape=guide_legend(override.aes = list(size = 3)))

#IL-6/JAK/STAT - ATC
library(qusage)
c2 <- read.gmt("c2.all.v2023.2.Hs.symbols.gmt")
atc_sig_all <- c(c2$MISIAK_ANAPLASTIC_THYROID_CARCINOMA_UP, c2$RODRIGUES_THYROID_CARCINOMA_ANAPLASTIC_UP) %>% unique()
jakstat_cebpb_regulon <- df_enriched %>% filter(Term == 'IL-6/JAK/STAT3 Signaling') %>% pull(Genes) %>% strsplit(., ", ")
epi <- AddModuleScore(epi, features = list(jakstat_cebpb_regulon[[1]]), name="jakstat_cebpb_regulon")
epi <- AddModuleScore(epi, features = list(atc_sig_all), name="atc_sig_all")

jakstat_atc <- data.frame(
  jakstat = epi$jakstat_cebpb_regulon1,
  atc = epi$atc_sig_all1,
  group = epi$annotation_epi)
jakstat_atc$group <- factor(jakstat_atc$group, levels = c('Follicular cell', 'PTC type1', 'PTC type2', 'ATC'))
jakstat_atc_cor <- cor.test(jakstat_atc$jakstat, jakstat_atc$atc)

ggplot(jakstat_atc, aes(x = jakstat, y = atc)) +
  geom_point(aes(shape = group, color = group)) +
  scale_color_manual(values = c('#EA8331', '#7CAE00', '#35A2FF', '#C00000')) +
  geom_smooth(method = 'lm', color ='red') +
  annotate("text", x = 0.9, y = 1.35, col = "black", size = 10,
           label = paste("r = ", signif(cor(jakstat_atc$jakstat, jakstat_atc$atc),3)))+
  theme_classic() +
  labs(x = 'JAK STAT Signature', y = 'ATC Signature') +
  theme(axis.title.x = element_text(size = 24, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 24, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 24, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 24, family = 'Arial', colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 24, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank()) &
  guides(shape=guide_legend(override.aes = list(size = 3)))


#Fig. 3D - IL-6/JAK/STAT, inflammatory response, and EMT gene set expression by epithelial cells in scRNA-seq data####
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

geneset_flt <- c(emt$gene,jakstat$gene,infla$gene)
gene_groups <- c(rep('Inflammatory response',7),rep('IL-6/JAK/STAT signaling',3),rep('EMT',9))

mat <- as.matrix(GetAssayData(epi, slot = "data")[geneset_flt, ])

mat <- t(scale(t(mat)))

mat[mat > 2.5] <- 2.5
mat[mat < -2.5] <- -2.5
meta <- epi@meta.data[colnames(mat), ]

set.seed(123)
cell_order <- unlist(lapply(split(rownames(meta), meta$annotation_epi), sample))

mat <- mat[, cell_order]
meta <- meta[cell_order, ]

anno_colors <- list(
  cell_type = setNames(c('#EA8331','#7CAE00', '#35A2FF', "#C00000"), unique(meta$annotation_epi))
)

ha <- HeatmapAnnotation(
  cell_type = meta$annotation_epi,
  col = anno_colors
)

col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),
  colors = c("navy", "white", "firebrick3")
)

Heatmap(
  mat,
  name = "Z-score",
  col = col_fun,
  split = gene_groups,
  top_annotation = ha,
  cluster_columns = FALSE,   # don’t cluster cells (optional)
  cluster_rows = TRUE,       # cluster genes
  show_column_names = FALSE, # hide cell names
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 23),   # gene name font size
  column_names_gp = gpar(fontsize = 23),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 23),  # legend title
    labels_gp = gpar(fontsize = 20)) ,
  row_title_gp = gpar(fontsize = 23) 
)


#Fig. 3E - IL-6/JAK/STAT, inflammatory response, and EMT gene set expression by Visium sample ####
ST_merge <- readRDS("R_objects/20251107_thca_ST_merge.rds")
ST_merge_tumor <- subset(ST_merge, celltype %in% c('Follicular cell','PTC type1','PTC type2','ATC'))
epi_DEG <- read.xlsx('Tumor_DEG/20231221_THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012/20240408_epi_nFeature_RNA250_v3_annotation_epi_DEG.xlsx', rowNames = TRUE)
infla <- epi_DEG %>% filter(cluster == 'ATC', avg_log2FC > 0.5, p_val_adj < 0.05, gene %in% infla_cebpb_regulon[[1]])
jakstat <- epi_DEG %>% filter(cluster == 'ATC', avg_log2FC > 0.5, p_val_adj < 0.05, gene %in% jakstat_cebpb_regulon[[1]])
emt <- epi_DEG %>% filter(cluster == 'ATC', avg_log2FC > 0.5, p_val_adj < 0.05, gene %in% emt_cebpb_regulon[[1]])

geneset_flt <- c(emt$gene,jakstat$gene,infla$gene)
gene_groups <- c(rep('EMT',9),rep('Inflammatory response',7),rep('JAK/STAT signaling',3))

expr <- FetchData(ST_merge_tumor, vars = c(geneset_flt, "sample"))
expr_scaled <- as.data.frame(scale(expr[, geneset_flt]))
expr_scaled$sample <- expr$sample
avg_scaled <- expr_scaled %>%
  group_by(sample) %>%
  summarise(across(all_of(geneset_flt), mean, na.rm = TRUE))

mat <- as.data.frame(avg_scaled)
rownames(mat) <- mat$sample
mat <- mat[, -1]  # remove sample column
mat <- t(as.matrix(mat))  # transpose to genes x samples
colnames(mat)[1:4] <- paste0('PT',1:4)
mat <- mat[c('EMP3','CD44','MMP14','ITGA5','HAS2','INHBA','ADM','NFKBIA','CSF2',
             'CD44','EMP3','TFPI2','MMP14','ITGA5','CTHRC1','INHBA',
             'CXCL3','VEGFA','LOX'),]
col_fun <- colorRamp2(c(-0.5, 0, 1), c("navy", "white", "firebrick3"))

group <- c(
  rep("PT", 4),     # PT1–PT4
  rep("PTC", 4),    # PTC1–PTC4
  rep("LPTC", 4),   # LPTC1–LPTC4
  rep("ATC", 4)             # ATC1 - ATC4
)
group <- factor(group, levels = c("PT", "PTC", "LPTC", "ATC"))

group_colors <- c(
  "PT" = "#4DAF4A",
  "PTC" = "#377EB8",
  "LPTC" = "#984EA3",
  "ATC" = "#E41A1C"
)

col_anno <- HeatmapAnnotation(
  CancerType = group,
  col = list(CancerType = group_colors),
  annotation_legend_param = list(
    CancerType = list(title = "Cancer type")
  )
)

Heatmap(
  mat,
  name = "Z-score",
  split = gene_groups,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = FALSE,
  top_annotation = col_anno,
  row_names_gp = gpar(fontsize = 23),   # gene name font size
  column_names_gp = gpar(fontsize = 23),
  heatmap_legend_param = list(
    title_position = "topcenter",
    title_gp = gpar(fontsize = 23),  # legend title
    labels_gp = gpar(fontsize = 18)),
  row_title_gp = gpar(fontsize = 23) 
  #cell_fun = function(j, i, x, y, width, height, fill) {
  #  if (mat[i, j] > 0) {
  #    grid.text("*", x, y, gp = gpar(fontsize = 12))
  #  }
  #}
)

#wilcox
mat_nodup <- mat[!duplicated(rownames(mat)), ]
sample <- colnames(mat_nodup)
group <- sub("[0-9]+$", "", sample)

atc_idx <- group == "ATC"

pvals <- apply(mat_nodup, 1, function(x) {
  wilcox.test(x[atc_idx], x[!atc_idx],
              alternative = "greater",
              exact = FALSE)$p.value
})

res <- data.frame(
  gene = rownames(mat_nodup),
  p_value = pvals,
  FDR = p.adjust(pvals, method = "BH"),
  row.names = NULL
)


#Fig. 3F - IL-6/JAK/STAT, EMT, inflammatory response signature in visium####
infla_cebpb_regulon <- df_enriched %>% filter(Term == 'Inflammatory Response') %>% pull(Genes) %>% strsplit(., ", ")
jakstat_cebpb_regulon <- df_enriched %>% filter(Term == 'IL-6/JAK/STAT3 Signaling') %>% pull(Genes) %>% strsplit(., ", ")
emt_cebpb_regulon <- df_enriched %>% filter(Term == 'Epithelial Mesenchymal Transition') %>% pull(Genes) %>% strsplit(., ", ")
DefaultAssay(ST_merge) <- 'Spatial'
ST_merge <- AddModuleScore(ST_merge, features = list(infla_cebpb_regulon[[1]]), name="infla_cebpb_regulon")
ST_merge <- AddModuleScore(ST_merge, features = list(jakstat_cebpb_regulon[[1]]), name="jakstat_cebpb_regulon")
ST_merge <- AddModuleScore(ST_merge, features = list(emt_cebpb_regulon[[1]]), name="emt_cebpb_regulon")

for (i in names(ST_merge@images)){
  assign(paste0('q_',i), SpatialFeaturePlot(ST_merge, c('emt_cebpb_regulon1','infla_cebpb_regulon1','jakstat_cebpb_regulon1'),
                                            pt.size.factor = 1.6, images = i, image.alpha = 0,  ncol = 3) &
           theme(text = element_text(family = "Arial")) & NoLegend())
}

plot_grid(q_N4,q_PTC2,q_LPTC4,q_ATC1, ncol = 1)


#Fig. 3G - IL-6/JAK/STAT, EMT, inflammatory response signature by visium samples####
library(ggplot2)
library(ggpubr)
library(tidyr)

expr <- FetchData(ST_merge, vars = c('infla_cebpb_regulon1','jakstat_cebpb_regulon1','emt_cebpb_regulon1', "sample", "group", "celltype"))
expr <- expr %>% filter(celltype %in% c('Follicular cell','PTC type1','PTC type2','ATC'))
expr$group <- factor(expr$group, levels = c("Normal", "PTC", "LPTC", "ATC"))
ATC_samples <- unique(expr$sample[expr$group == "ATC"])
signatures <- c("infla_cebpb_regulon1", "jakstat_cebpb_regulon1", "emt_cebpb_regulon1")
test_results <- lapply(signatures, function(reg) {
  tmp <- lapply(ATC_samples, function(s) {
    atc_expr <- expr %>% filter(sample == s)
    others_expr <- expr %>% filter(group != "ATC")
    
    test <- wilcox.test(atc_expr[[reg]], others_expr[[reg]])
    
    data.frame(
      regulon = reg,
      group1 = "Others",
      group2 = s,
      p_value = test$p.value
    )
  }) %>% bind_rows()
  tmp
}) %>% bind_rows()

test_results <- test_results %>%
  group_by(regulon) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

test_results$label <- ifelse(test_results$p_adj < 0.001, "***",
                             ifelse(test_results$p_adj < 0.01, "**",
                                    ifelse(test_results$p_adj < 0.05, "*", "ns")))

expr_long <- expr %>%
  filter(celltype %in% c('Follicular cell','PTC type1','PTC type2','ATC')) %>%
  select(infla_cebpb_regulon1, jakstat_cebpb_regulon1, emt_cebpb_regulon1, sample, group) %>%
  pivot_longer(
    cols = c(infla_cebpb_regulon1, jakstat_cebpb_regulon1, emt_cebpb_regulon1),
    names_to = "regulon",
    values_to = "value"
  )

expr_long$sample <- gsub('^N','PT',expr_long$sample)
expr_long$sample <- factor(expr_long$sample, levels = c(paste0('PT',1:4),paste0('PTC',1:4),paste0('LPTC',1:4),paste0('ATC',1:4)))
expr_long$group <- gsub('^Normal','PT',expr_long$group)
expr_long$group <- factor(expr_long$group, levels = c("PT", "PTC", "LPTC", "ATC"))

p_emt <- ggplot(expr_long[expr_long$regulon == 'emt_cebpb_regulon1',], aes(x = sample, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.5) +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.25) +
  scale_fill_manual(values = c("PT" = "#4DAF4A",
                               "PTC" = "#377EB8",
                               "LPTC" = "#984EA3",
                               "ATC" = "#E41A1C")) +
  theme_classic() +
  NoLegend() +
  labs(title = "EMT\nsignature score") +
  theme(title = element_text(size = 21, color = 'black'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 21,angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 21,color = 'black'))

p_inflam <- ggplot(expr_long[expr_long$regulon == 'infla_cebpb_regulon1',], aes(x = sample, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.5) +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.25) +
  scale_fill_manual(values = c("PT" = "#4DAF4A",
                               "PTC" = "#377EB8",
                               "LPTC" = "#984EA3",
                               "ATC" = "#E41A1C")) +
  theme_classic() +
  NoLegend() +
  labs(title = "Inflammatory response\nsignature score") +
  theme(title = element_text(size = 21, color = 'black'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 21,angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 21,color = 'black'))

p_jakstat <- ggplot(expr_long[expr_long$regulon == 'jakstat_cebpb_regulon1',], aes(x = sample, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.5) +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.25) +
  scale_fill_manual(values = c("PT" = "#4DAF4A",
                               "PTC" = "#377EB8",
                               "LPTC" = "#984EA3",
                               "ATC" = "#E41A1C")) +
  theme_classic() +
  NoLegend() +
  labs(title = "IL-6/JAK/STAT\nsignature score") +
  theme(title = element_text(size = 21, color = 'black'),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 21,angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 21,color = 'black'))

plot_grid(p_emt,p_inflam,p_jakstat, ncol = 3)


#Fig. 3H - pathways enriched in CEBPB overexpressed PTC cell line####
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)

counts <- as.matrix(read.csv("bulk/count_matrix/gene_count_matrix.csv",sep=",",row.names="gene_id"))
counts <- counts[,c(4,5,6,1,2,3)]

sample <- colnames(counts)
condition <- factor(rep(c("CON","CEBPB"),c(3,3) ))
metadata <- data.frame(row.names = sample,condition)
all(colnames(counts) %in% rownames(metadata))
all(rownames(metadata) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds$condition <- relevel(dds$condition, ref = "CON")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "CEBPB", "CON"))

resLFC <- lfcShrink(dds, coef="condition_CEBPB_vs_CON" , type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
resLFC <- as.data.frame(resLFC)
row.names(resLFC) <- gsub('\\.[0-9]*$', '', row.names(resLFC))
resLFC$Ensembl<-row.names(resLFC)

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=resLFC$Ensembl, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
resLFC <- inner_join(resLFC, ens2symbol, by=c("Ensembl"="ENSEMBL"))

resLFC_gsea <- resLFC %>% 
  dplyr::select(SYMBOL, log2FoldChange) %>%
  na.omit() %>% 
  distinct() %>%
  group_by(SYMBOL) %>% 
  summarize(stat=mean(log2FoldChange))

ranks <- deframe(resLFC_gsea)
set.seed(0)
fgseaRes <- fgsea(pathways=fgsea_hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy_flt <- fgseaResTidy %>% filter(padj < 0.05)
fgseaResTidy_flt$pathway <- gsub('HALLMARK_','',fgseaResTidy_flt$pathway)
colnames(fgseaResTidy_flt)[3] <- 'Adjusted P.value'

ggplot(fgseaResTidy_flt, 
       aes(x = NES, y = reorder(pathway, NES), fill = `Adjusted P.value`)) +
  geom_bar(stat = 'identity') +
  scale_fill_viridis(direction = -1) +
  theme_classic() +
  labs(y = "", x = "NES") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18, color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "top", 
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, color = "black"),
        axis.ticks = element_line(color = "black", size=2))
