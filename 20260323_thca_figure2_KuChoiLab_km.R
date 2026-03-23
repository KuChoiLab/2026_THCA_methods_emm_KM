setwd('/Users/yoo/Dropbox/THCA_scrna/')
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

# Fig. 2A - Epithelial cell UMAP ####
epi_sub3 <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')
DimPlot(epi_sub3, group.by = 'annotation_epi', cols = c('#EA8331','#7CAE00', '#35A2FF', "#C00000")) & 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15, family="Arial"),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) &
  guides(col=guide_legend(ncol=1, override.aes = list(size = 3))) &
  NoAxes() &
  ggtitle(NULL)


# Fig. 2B - TDS ####
epi_sub3 <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')
TDS <- c('TG', 'TPO', 'SLC26A4', 'DIO2', 'TSHR', 'PAX8', 'DUOX1', 'DUOX2', 'NKX2-1', 'GLIS3', 'FOXE1', 'TFF3', 'FHL1')
epi_sub3 <- AddModuleScore(epi_sub3, features = list(TDS), name="TDS")

df <- epi_sub3@meta.data

pairs <- list(
  c("Follicular cell", "PTC type1"),
  c("Follicular cell", "PTC type2"),
  c("ATC", "PTC type1"),
  c("ATC", "PTC type2")
)

pvals <- sapply(pairs, \(p)
                signif(wilcox.test(df$TDS1[df$annotation_epi_v2 == p[1]],
                                   df$TDS1[df$annotation_epi_v2 == p[2]])$p.value, 4)
)

library(ggplot2)
library(ggsignif)

ggplot(epi_sub3@meta.data, aes(x=annotation_epi_v2, y=TDS1, fill=annotation_epi_v2)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() + 
  NoLegend() +
  labs(x = 'Cell type', y = 'TDS') +
  theme(title = element_blank(),
        axis.title.x = element_text(size = 0, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 28, family = 'Arial', colour = 'black', angle = 30, hjust=1),
        axis.title.y = element_text(size = 28, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 28, family = 'Arial', colour = 'black')) &
  scale_fill_manual(values = c('#EA8331','#7CAE00', '#35A2FF', "#C00000")) &
  geom_signif(
    y_position = c(2,1.8,1.6,1.4), xmin = c(1,1,2,3), xmax = c(2,3,4,4),
    annotation = c('***', '***', '***', '***'), 
    tip_length = 0.02, textsize = 10, color = 'black'
  )


# Fig. 2C - epithelial cell marekers ####
epi <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')

DotPlot(epi, features = c('TG','TPO','TFF3','FN1','GDF15','SLC34A2', 'S100A1', 'EEF1A1','RPS4X', 'KRT19', 'UBE2C', 'S100A10', 'TFPI2'), 
        group.by = 'annotation_epi', cols = c('white', 'red')) +
  RotatedAxis() +
  coord_trans() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(colour = "black", angle = 30, size = 20, vjust = 0.9, hjust = 0.95),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(),
        panel.border = element_blank())


# Fig. 2D - epithelial cell marekers ####
##nebula####
epi <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')
epi@meta.data <- epi@meta.data %>% mutate(annotation_epi_v2 = case_when(
  annotation_epi == 'Follicular cell' ~ 'Follicular cells',
  annotation_epi == 'PTC type1' ~ 'PTC type1',
  annotation_epi == 'PTC type2' ~ 'PTC type2',
  annotation_epi == 'ATC' ~ 'ATC',
  TRUE ~ 'Others'
))
epi@meta.data$annotation_epi_v2 <- factor(epi@meta.data$annotation_epi_v2, levels = c('Follicular cells', 'PTC type1', 'PTC type2', 'ATC'))

for (i in unique(epi$annotation_epi_v2)) {
  epi@meta.data <- epi@meta.data %>% dplyr::mutate(nebula_anno = case_when(
    annotation_epi_v2 == i ~ paste0(i),
    TRUE ~ "Others"))
  epi$nebula_anno <- factor(epi$nebula_anno, levels = rev(c(paste0(i), "Others"))) # 중요
  seuratdata <- scToNeb(obj = epi, assay = "RNA", id = "orig.ident", pred = c("nebula_anno", "orig.ident"), offset = "nCount_RNA")
  df = model.matrix(~nebula_anno+orig.ident, data = seuratdata$pred)
  re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)",paste0("nebula_anno",i))],offset=seuratdata$offset)
  saveRDS(re, paste0('nebula/epi/nebula_',gsub(' ','_',i),'_km.rds'))
  res <- re$summary %>% dplyr::select(gene, paste0('logFC_nebula_anno',i), paste0('se_nebula_anno',i),paste0('p_nebula_anno',i))
  res$Padj <- p.adjust(res[,paste0('p_nebula_anno',i)], method = "BH")
  colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
  write.xlsx(res,paste0('nebula/epi/20260127_thca_scrna_nebula_',gsub(' ','_',i),'_km.xlsx'))}

res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(20)

res_all_top100 <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("gene", "avg_log2FC", "se", "p_val", "p_val_adj"))
for (i in unique(epi$annotation_epi_v2)) {
  res <- read.xlsx(paste0('nebula/epi/20260127_thca_scrna_nebula_',gsub(' ','_',i),'_km.xlsx'))
  res$p_val_adj <- p.adjust(res[,'p_val'], method = "BH")
  res$celltype <- paste0(i)
  res <- res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(100)
  res_all_top100 <- rbind(res_all_top100,res)}
write.xlsx(res_all_top100,"nebula/epi/20260127_thca_scrna_nebula_Epi_top100DEGs_km.xlsx")

res_all_top100 %>% filter(celltype == 'ATC', p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(20)


## GSEA ####
library(fgsea)
library(tibble)

# Follicular cells
celltype <- "Follicular cells"
res <- read.xlsx(paste0('nebula/epi/20260127_thca_scrna_nebula_',gsub(' ','_',celltype),'_km.xlsx'))
res$p_val_adj <- p.adjust(res[,'p_val'], method = "BH")
df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)
df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
set.seed(1234); fgseaRes <- fgseaMultilevel(pathways = pw, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- as.data.frame(fgseaRes)
fgseaRes_final_df <- data.frame()
for (i in seq(1:nrow(fgseaRes))) {
  fgseaRes_final_df_tmp <- fgseaRes[i,]
  leading_edges_tmp <- fgseaRes[i,]$leadingEdge %>% unlist()
  fgseaRes_final_df_tmp$leadingEdge <- paste(leading_edges_tmp, collapse = "/")
  
  fgseaRes_final_df <- rbind(fgseaRes_final_df, fgseaRes_final_df_tmp)
}
fgseaRes_final_df$celltype <- celltype
write.xlsx(fgseaRes_final_df,"nebula/epi/20260127_thca_scrna_nebula_Follicularcells_GSEA_GOBP_km.xlsx")
fgseaRes_final_df %>% filter(padj < 0.05) %>% arrange(desc(NES)) %>% select(pathway,NES) %>% head(20)


# PTC1
celltype <- "PTC type1"
res <- read.xlsx(paste0('nebula/epi/20260127_thca_scrna_nebula_',gsub(' ','_',celltype),'_km.xlsx'))
res$p_val_adj <- p.adjust(res[,'p_val'], method = "BH")
df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)
df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
set.seed(1234); fgseaRes <- fgseaMultilevel(pathways = pw, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- as.data.frame(fgseaRes)
fgseaRes_final_df <- data.frame()
for (i in seq(1:nrow(fgseaRes))) {
  fgseaRes_final_df_tmp <- fgseaRes[i,]
  leading_edges_tmp <- fgseaRes[i,]$leadingEdge %>% unlist()
  fgseaRes_final_df_tmp$leadingEdge <- paste(leading_edges_tmp, collapse = "/")
  
  fgseaRes_final_df <- rbind(fgseaRes_final_df, fgseaRes_final_df_tmp)
}
fgseaRes_final_df$celltype <- celltype
write.xlsx(fgseaRes_final_df,"nebula/epi/20260127_thca_scrna_nebula_PTCtype1_GSEA_GOBP_km.xlsx")

# PTC2
celltype <- "PTC type2"
res <- read.xlsx(paste0('nebula/epi/20260127_thca_scrna_nebula_',gsub(' ','_',celltype),'_km.xlsx'))
res$p_val_adj <- p.adjust(res[,'p_val'], method = "BH")
df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)
df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
set.seed(1234); fgseaRes <- fgseaMultilevel(pathways = pw, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- as.data.frame(fgseaRes)
fgseaRes_final_df <- data.frame()
for (i in seq(1:nrow(fgseaRes))) {
  fgseaRes_final_df_tmp <- fgseaRes[i,]
  leading_edges_tmp <- fgseaRes[i,]$leadingEdge %>% unlist()
  fgseaRes_final_df_tmp$leadingEdge <- paste(leading_edges_tmp, collapse = "/")
  
  fgseaRes_final_df <- rbind(fgseaRes_final_df, fgseaRes_final_df_tmp)
}
fgseaRes_final_df$celltype <- celltype
write.xlsx(fgseaRes_final_df,"nebula/epi/20260127_thca_scrna_nebula_PTCtype2_GSEA_GOBP_km.xlsx")
fgseaRes_final_df %>% filter(padj < 0.05) %>% arrange(desc(NES)) %>% select(pathway,NES) %>% head(20)


# ATC
celltype <- "ATC"
res <- read.xlsx(paste0('nebula/epi/20260127_thca_scrna_nebula_',gsub(' ','_',celltype),'_km.xlsx'))
res$p_val_adj <- p.adjust(res[,'p_val'], method = "BH")
df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
set.seed(1234); fgseaRes <- fgseaMultilevel(pathways = pw, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- as.data.frame(fgseaRes)
fgseaRes_final_df <- data.frame()
for (i in seq(1:nrow(fgseaRes))) {
  fgseaRes_final_df_tmp <- fgseaRes[i,]
  leading_edges_tmp <- fgseaRes[i,]$leadingEdge %>% unlist()
  fgseaRes_final_df_tmp$leadingEdge <- paste(leading_edges_tmp, collapse = "/")
  
  fgseaRes_final_df <- rbind(fgseaRes_final_df, fgseaRes_final_df_tmp)
}
fgseaRes_final_df$celltype <- celltype
write.xlsx(fgseaRes_final_df,"nebula/epi/20260127_thca_scrna_nebula_ATC_GSEA_GOBP_km.xlsx")


#Plot
epi_sub3_follicular_DEG_gobp <- read.xlsx('nebula/epi/20260127_thca_scrna_nebula_Follicularcells_GSEA_GOBP_km.xlsx')
epi_sub3_PTCtype1_DEG_gobp <- read.xlsx('nebula/epi/20260127_thca_scrna_nebula_PTCtype1_GSEA_GOBP_km.xlsx')
epi_sub3_PTCtype2_DEG_gobp <- read.xlsx('nebula/epi/20260127_thca_scrna_nebula_PTCtype2_GSEA_GOBP_km.xlsx')
epi_sub3_ATC_DEG_gobp <- read.xlsx('nebula/epi/20260127_thca_scrna_nebula_ATC_GSEA_GOBP_km.xlsx')

epi_sub3_folli_DEG_gobp_rep_terms <- c('thyroid hormone metabolic process','modified amino acid metabolic process', 'antigen processing and presentation of peptide antigen','oxidative phosphorylation','cellular response to metal ion')
epi_sub3_PTCtype1_DEG_gobp_rep_terms <- c('response to growth hormone','cell junction assembly', 'regulation of mrna metabolic process',
                                          'positive regulation of cell substrate adhesion','stem cell differentiation')
epi_sub3_PTCtype2_DEG_gobp_rep_terms <- c('ribosome biogenesis','cytoplasmic translation', 'myeloid leukocyte activation',
                                          'rRNA metabolic process','immune effector process')
epi_sub3_ATC_DEG_gobp_rep_terms <- c('positive regulation of cell cycle checkpoint','regulation of cytoskeleton organization', 'DNA replication',
                                     'regulation of cell-cell adhesion','activation of immune response')

enrichmate_all_stat_rep <- c(epi_sub3_folli_DEG_gobp_rep_terms,epi_sub3_PTCtype1_DEG_gobp_rep_terms,epi_sub3_PTCtype2_DEG_gobp_rep_terms,epi_sub3_ATC_DEG_gobp_rep_terms)

enrichmate_all_stat_rep <- gsub('^','GOBP_',enrichmate_all_stat_rep)
enrichmate_all_stat_rep <- gsub('-','_',enrichmate_all_stat_rep)
enrichmate_all_stat_rep <- gsub(' ','_',enrichmate_all_stat_rep)
enrichmate_all_stat_rep <- toupper(enrichmate_all_stat_rep)
enrichmate_all_stat_rep <- gsub(',','',enrichmate_all_stat_rep)

epi_sub3_follicular_DEG_gobp_flt <- epi_sub3_follicular_DEG_gobp %>% filter(pathway %in% enrichmate_all_stat_rep)
epi_sub3_PTCtype1_DEG_gobp_flt <- epi_sub3_PTCtype1_DEG_gobp %>% filter(pathway %in% enrichmate_all_stat_rep)
epi_sub3_PTCtype2_DEG_gobp_flt <- epi_sub3_PTCtype2_DEG_gobp %>% filter(pathway %in% enrichmate_all_stat_rep)
epi_sub3_ATC_DEG_gobp_flt <- epi_sub3_ATC_DEG_gobp %>% filter(pathway %in% enrichmate_all_stat_rep)

epi_sub3_follicular_DEG_gobp_flt$group <- 'Follicular cells'
epi_sub3_PTCtype1_DEG_gobp_flt$group <- 'PTC type1'
epi_sub3_PTCtype2_DEG_gobp_flt$group <- 'PTC type2'
epi_sub3_ATC_DEG_gobp_flt$group <- 'ATC'

epi_sub3_all_DEG_gobp_flt <- do.call(rbind, list(epi_sub3_follicular_DEG_gobp_flt,epi_sub3_PTCtype1_DEG_gobp_flt,epi_sub3_PTCtype2_DEG_gobp_flt,epi_sub3_ATC_DEG_gobp_flt))
epi_sub3_all_DEG_gobp_flt2 <- epi_sub3_all_DEG_gobp_flt %>% filter(padj < 0.05)
epi_sub3_all_DEG_gobp_flt2$group <- factor(epi_sub3_all_DEG_gobp_flt2$group, levels = c('Follicular cells', 'PTC type1', 'PTC type2', 'ATC'))
epi_sub3_all_DEG_gobp_flt2$pathway <- gsub('^GOBP_','',epi_sub3_all_DEG_gobp_flt2$pathway)
epi_sub3_all_DEG_gobp_flt2$pathway <- factor(epi_sub3_all_DEG_gobp_flt2$pathway, levels = gsub('^GOBP_','',enrichmate_all_stat_rep))

ggplot(epi_sub3_all_DEG_gobp_flt2, aes(group,pathway)) + geom_point(aes(color = NES, size = -log10(padj))) + 
  ylab(NULL) + theme_minimal() +
  theme(axis.text.x = element_text(colour = "black", angle = 60, size = 18, vjust=1, hjust=1), 
        axis.text.y = element_text(colour = "black", size = 18), plot.title = element_blank(), 
        legend.position = "right", legend.text = element_text(size = 18),
        legend.title = element_text(size = 18), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.spacing = unit(1, 'cm'), legend.spacing.x = unit(0.5,'cm')) + 
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) + 
  scale_size(range = c(1, 12), breaks = c(0, 3, 6, 9, 12)) + scale_color_gradientn(colours=c("blue","pink","red")) +
  scale_shape_manual(values = c(13,16))


#Fig. 2E - infercnv####
library(Seurat)
library(dplyr)
library(infercnv)
setwd('/Users/yoo/Dropbox/THCA_scrna/')

epi_mnn <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')
myeloid <- readRDS('R_objects/THCA_myeloid_mnn_20231019.rds')
myeloid_sub <- subset(myeloid, myeloid_sub %in% c('MT contam (myeloid)', 'T contam (myeloid)'), invert =T)
myeloid_sub$myeloid_sub %>% table
Mye_thy <- merge(myeloid_sub, epi_mnn)

Mye_thy@meta.data <- Mye_thy@meta.data %>% mutate(infercnv_group = case_when(
  annotation_epi == "Follicular cell" ~ "Follicular cell", 
  annotation_epi == "PTC type1" ~ "PTC type1",
  annotation_epi == "PTC type2" ~ "PTC type2",
  annotation_epi == "ATC" ~ "ATC",
  myeloid_sub %in% c("C1QC+ TAM", "CD16+ Monocytes", "cDC1", "cDC2", "cDC3", "Classical Monocytes", "IL-8+ Monocyte", "Mast cells", "pDC",
                     "Proliferative Myeloid", "SPP1+ TAM") ~ "Myeloid cell", 
  TRUE ~ "Others"))
Mye_thy_counts <- as.matrix(Mye_thy@assays$RNA@counts)

#annotations_file
Mye_thy@meta.data %>% select(infercnv_group) -> Mye_thy_Annotations
Mye_thy_Annotations[,2] <- row.names(Mye_thy_Annotations)
Mye_thy_Annotations[,c("V2", "infercnv_group")] -> Mye_thy_Annotations
dir.create('infercnv_epi_sub3')
write.table(Mye_thy_Annotations, file = "/Users/yoo/Dropbox/THCA_scrna/infercnv_epi_sub3/Mye_thy_Annotations.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#gene_order_pos
library(Seurat)
library(SeuratDisk)
SaveH5Seurat(Mye_thy, filename = "infercnv_epi_sub3/Mye_thy.h5Seurat")
Convert("infercnv_epi_sub3/Mye_thy.h5Seurat", dest = "h5ad")
#after generating Mye_thy_pos.csv in jupyter notebook (macbook conda infercnvpy, 20240414_thca_epi_sub3.ipynb)
Mye_thy_pos <- read.table("infercnv_epi_sub3/Mye_thy_pos.csv", sep = "\t")
Mye_thy_pos <- Mye_thy_pos[-1,]
Mye_thy_pos_v2 <- na.omit(Mye_thy_pos)
Mye_thy_pos_final <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("V1","V2", "V3", "V4"))
chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
         "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
for (i in chr) {
  split <- Mye_thy_pos_v2 %>% filter(V2 == i) %>% arrange(as.numeric(V3))
  Mye_thy_pos_final <- rbind(Mye_thy_pos_final, split)
}
head(Mye_thy_pos_final)
tail(Mye_thy_pos_final)
write.table(Mye_thy_pos_final, file = "infercnv_epi_sub3/Mye_thy_pos.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

Mye_thy_infercnv_obj = CreateInfercnvObject(raw_counts_matrix=Mye_thy_counts,
                                            annotations_file="infercnv_epi_sub3/Mye_thy_Annotations.txt",
                                            delim="\t",
                                            gene_order_file="infercnv_epi_sub3/Mye_thy_pos.txt",
                                            ref_group_names=c("Myeloid cell"))

Mye_thy_infercnv_obj = infercnv::run(Mye_thy_infercnv_obj,
                                     cutoff=0.1,
                                     out_dir="~/Dropbox/THCA_scrna/infercnv_epi_sub3/out_Mye_thyrocyte_PTC_ATC",
                                     cluster_by_groups=T,
                                     plot_steps=F,
                                     denoise=T,
                                     sd_amplifier=3,
                                     noise_logistic=TRUE,
                                     HMM=F)

##re-run with highlighted Myeloid
infercnvR_Myeloid_thyrocyte_references <- read.table("infercnv_epi_sub3/out_Mye_thyrocyte_PTC_ATC/infercnv.references.txt",sep=" ",header=TRUE)
Mye_thy_infercnv_obj <- readRDS('infercnv_epi_sub3/out_Mye_thyrocyte_PTC_ATC/run.final.infercnv_obj')
chr6_genes <- Mye_thy_infercnv_obj@gene_order %>% filter(chr == "chr6")
infercnvR_Myeloid_thyrocyte_references_chr6 <- infercnvR_Myeloid_thyrocyte_references[row.names(infercnvR_Myeloid_thyrocyte_references) %in% row.names(chr6_genes),]
infercnvR_Myeloid_thyrocyte_references_chr6_mean <- apply(infercnvR_Myeloid_thyrocyte_references_chr6, 2, mean) #colmean
length(which(infercnvR_Myeloid_thyrocyte_references_chr6_mean < 1)) 
length(which(infercnvR_Myeloid_thyrocyte_references_chr6_mean > 1)) 
row.names(infercnvR_Myeloid_thyrocyte_references_chr6)[114:129] #114-129 -> MHC
apply(infercnvR_Myeloid_thyrocyte_references_chr6[114:129,5000:5023],2,mean) #identifying del region 
highlighted_mye <- colnames(infercnvR_Myeloid_thyrocyte_references_chr6)[3033:5023]
highlighted_mye_v2 <- gsub(".", "-", highlighted_mye, fixed = TRUE)
#subsetting highlighted cells raw_counts_matrix
epi_mnn <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')
myeloid <- readRDS('R_objects/THCA_myeloid_mnn_20231019.rds')
myeloid_sub <- subset(myeloid, myeloid_sub %in% c('MT contam (myeloid)', 'T contam (myeloid)'), invert =T)

myeloid_highlight <- subset(myeloid_sub, cells= highlighted_mye_v2)
thy_Mye_highlighted <- merge(myeloid_highlight, epi_mnn)
thy_Mye_highlighted@meta.data <- thy_Mye_highlighted@meta.data %>% mutate(infercnv_group = case_when(
  annotation_epi == "Follicular cell" ~ "Follicular cell", 
  annotation_epi == "PTC type1" ~ "PTC type1",
  annotation_epi == "PTC type2" ~ "PTC type2",
  annotation_epi == "ATC" ~ "ATC",
  myeloid_sub %in% c("C1QC+ TAM", "CD16+ Monocytes", "cDC1", "cDC2", "cDC3", "Classical Monocytes", "IL-8+ Monocyte", "Mast cells", "pDC",
                     "Proliferative Myeloid", "SPP1+ TAM") ~ "Myeloid cell", 
  TRUE ~ "Others"))

thy_Mye_highlighted_counts <- as.matrix(thy_Mye_highlighted@assays$RNA@counts)

#annotations_file
thy_Mye_highlighted@meta.data %>% select(infercnv_group) -> thy_Mye_highlighted_Annotations
thy_Mye_highlighted_Annotations[,2] <- row.names(thy_Mye_highlighted_Annotations)
thy_Mye_highlighted_Annotations[,c("V2", "infercnv_group")] -> thy_Mye_highlighted_Annotations
write.table(thy_Mye_highlighted_Annotations, file = "/Users/yoo/Dropbox/THCA_scrna/infercnv_epi_sub3/thy_Mye_highlighted_Annotations.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#gene_order_pos
library(Seurat)
library(SeuratDisk)
SaveH5Seurat(thy_Mye_highlighted, filename = "infercnv_epi_sub3/thy_Mye_highlighted.h5Seurat")
Convert("infercnv_epi_sub3/thy_Mye_highlighted.h5Seurat", dest = "h5ad")
#after generating thy_Mye_highlighted_pos.csv in jupyter notebook (macbook conda infercnvpy, 20240414_thca_epi_sub3_myeloid_highlighted)
thy_Mye_highlighted_pos <- read.table("infercnv_epi_sub3/Mye_thy_highlighted_pos.csv", sep = "\t")
thy_Mye_highlighted_pos <- thy_Mye_highlighted_pos[-1,]
thy_Mye_highlighted_pos_v2 <- na.omit(thy_Mye_highlighted_pos)
thy_Mye_highlighted_pos_final <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("V1","V2", "V3", "V4"))
chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
         "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
for (i in chr) {
  split <- thy_Mye_highlighted_pos_v2 %>% filter(V2 == i) %>% arrange(as.numeric(V3))
  thy_Mye_highlighted_pos_final <- rbind(thy_Mye_highlighted_pos_final, split)
}
write.table(thy_Mye_highlighted_pos_final, file = "infercnv_epi_sub3/thy_Mye_highlighted_pos.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#run infercnv in imac
library(infercnv)
thy_Mye_highlighted_infercnv_obj = CreateInfercnvObject(raw_counts_matrix=thy_Mye_highlighted_counts,
                                                        annotations_file="infercnv_epi_sub3/thy_Mye_highlighted_Annotations.txt",
                                                        delim="\t",
                                                        gene_order_file="infercnv_epi_sub3/thy_Mye_highlighted_pos.txt",
                                                        ref_group_names=c("Myeloid cell"))

thy_Mye_highlighted_infercnv_obj = infercnv::run(thy_Mye_highlighted_infercnv_obj,
                                                 cutoff=0.1,
                                                 out_dir="~/Dropbox/THCA_scrna/infercnv_epi_sub3/out_thy_Mye_PTC_ATC_highlighted",
                                                 cluster_by_groups=T,
                                                 plot_steps=F,
                                                 denoise=T,
                                                 sd_amplifier=3,
                                                 noise_logistic=TRUE,
                                                 HMM=F)


#Fig. 2F - CEBPB featureplot####
epi_sub3 <- readRDS('R_objects/THCA_epithelial_no_PTCEO_mnn_flt_mnn_20231012_nFeature_RNA250_v3.rds')
FeaturePlot(epi_sub3, 'CEBPB', cols = c('grey', 'red')) & 
  theme(legend.position = "right",
        legend.text = element_text(size = 15, family="Arial"),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) &
  guides(fill=guide_legend(ncol=1, override.aes = list(size = 3))) &
  NoAxes() &
  ggtitle(NULL)


#Fig. 2G - CEBPB expression in Visium (boxplot)####
ST_merge <- readRDS("R_objects/20251107_thca_ST_merge.rds")
ST_merge_tumor <- subset(ST_merge, celltype %in% c('Follicular cell','PTC type1','PTC type2','ATC'))

cebpb_expression <- FetchData(ST_merge_tumor, vars = "CEBPB")
cebpb_expression$group <- ST_merge_tumor$group
cebpb_expression$sample <- ST_merge_tumor$sample
cebpb_expression$group <- gsub('Normal','PT',cebpb_expression$group)
cebpb_expression$group <- factor(cebpb_expression$group, levels = c('PT','PTC','LPTC','ATC'))
cebpb_expression$sample <- gsub('N','PT',cebpb_expression$sample)
cebpb_expression$sample <- factor(cebpb_expression$sample, levels = c(paste0('PT',1:4),paste0('PTC',1:4),paste0('LPTC',1:4),paste0('ATC',1:4)))

ggplot(cebpb_expression, aes(x = sample, y = CEBPB, fill = group)) +
  geom_violin() +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  ylim(0,5) +
  theme_classic() +
  NoLegend() +
  #labs(y = 'CEBP/β') +
  theme(axis.text.x = element_text(colour = "black", angle = 60, size = 22, vjust=1, hjust=1), 
        axis.text.y = element_text(colour = "black", size = 22), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(colour = "black", size = 22),
        #plot.title = element_text(colour = "black", size = 15, hjust = 0.5), 
        #legend.position = "right", legend.text = element_text(size = 15),
        #legend.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size=0.5),
        #legend.spacing = unit(1, 'cm'), legend.spacing.x = unit(0.5,'cm')
  ) +
  scale_fill_manual(values = c("PT" = "#4DAF4A",
                               "PTC" = "#377EB8",
                               "LPTC" = "#984EA3",
                               "ATC" = "#E41A1C"))


#nebula
library(nebula)
ST_merge$group <- gsub('Normal','PT',ST_merge$group)
ST_merge$sample <- gsub('^N','PT',ST_merge$sample)
ST_merge$sample <- factor(ST_merge$sample, levels = c(paste0('PT',1:4),paste0('PTC',1:4),paste0('LPTC',1:4),paste0('ATC',1:4)))
ST_merge_tumor <- subset(ST_merge, celltype %in% c('Follicular cell','PTC type1','PTC type2','ATC'))
##PT vs. ATC
ST_merge_atc_pt <- subset(ST_merge_tumor, group %in% c('ATC','PT'))
ST_merge_atc_pt$group <- factor(ST_merge_atc_pt$group, levels = rev(c("ATC", "PT"))) # 중요
ST_merge_atc_pt$sample <- droplevels(ST_merge_atc_pt$sample)
seuratdata <- scToNeb(obj = ST_merge_atc_pt, assay = "Spatial", id = "sample", pred = c("group", "sample"), offset = "nCount_Spatial")
df = model.matrix(~group+sample, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
#data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 CTL")],offset=seuratdata$offset)
#re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)
re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","groupATC")],offset=seuratdata$offset)

saveRDS(re, "visium_public/nebula/nebula_atc_pt_km.rds")
re <- readRDS("visium_public/nebula/nebula_atc_pt_km.rds")
res <- re$summary %>% dplyr::select(gene, logFC_groupATC, se_groupATC, p_groupATC)
res$Padj <- p.adjust(res$p_groupATC, method = "BH")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
res[res$gene == 'CEBPB',]
write.xlsx(res,"visium_public/nebula/20260127_thca_scrna_nebula_atc_pt_km.xlsx")

##PTC vs. ATC
ST_merge_atc_ptc <- subset(ST_merge_tumor, group %in% c('ATC','PTC'))
ST_merge_atc_ptc$group <- factor(ST_merge_atc_ptc$group, levels = rev(c("ATC", "PTC"))) # 중요
ST_merge_atc_ptc$sample <- droplevels(ST_merge_atc_ptc$sample)
seuratdata <- scToNeb(obj = ST_merge_atc_ptc, assay = "Spatial", id = "sample", pred = c("group", "sample"), offset = "nCount_Spatial")
df = model.matrix(~group+sample, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
#data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 CTL")],offset=seuratdata$offset)
#re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)
re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","groupATC")],offset=seuratdata$offset)

saveRDS(re, "visium_public/nebula/nebula_atc_ptc_km.rds")
re <- readRDS("visium_public/nebula/nebula_atc_ptc_km.rds")
res <- re$summary %>% dplyr::select(gene, logFC_groupATC, se_groupATC, p_groupATC)
res$Padj <- p.adjust(res$p_groupATC, method = "BH")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
res[res$gene == 'CEBPB',]
write.xlsx(res,"visium_public/nebula/20260127_thca_scrna_nebula_atc_ptc_km.xlsx")

##LPTC vs. ATC
ST_merge_atc_lptc <- subset(ST_merge_tumor, group %in% c('ATC','LPTC'))
ST_merge_atc_lptc$group <- factor(ST_merge_atc_lptc$group, levels = rev(c("ATC", "LPTC"))) # 중요
ST_merge_atc_lptc$sample <- droplevels(ST_merge_atc_lptc$sample)
seuratdata <- scToNeb(obj = ST_merge_atc_lptc, assay = "Spatial", id = "sample", pred = c("group", "sample"), offset = "nCount_Spatial")
df = model.matrix(~group+sample, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
#data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 CTL")],offset=seuratdata$offset)
#re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)
re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","groupATC")],offset=seuratdata$offset)

saveRDS(re, "visium_public/nebula/nebula_atc_lptc_km.rds")
re <- readRDS("visium_public/nebula/nebula_atc_lptc_km.rds")
res <- re$summary %>% dplyr::select(gene, logFC_groupATC, se_groupATC, p_groupATC)
res$Padj <- p.adjust(res$p_groupATC, method = "BH")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
res[res$gene == 'CEBPB',]
write.xlsx(res,"visium_public/nebula/20260127_thca_scrna_nebula_atc_lptc_km.xlsx")


#Fig. 2H - CEBPB expression in Visium (Featureplot)####
ST_merge <- readRDS("R_objects/20251107_thca_ST_merge.rds")

SpatialFeaturePlot(
  ST_merge,
  features = "CEBPB",
  images = c('N2','PTC2','LPTC3','ATC1'),
  ncol = 4
) & scale_fill_viridis_c(option = "magma") & NoLegend()


#Fig. 2I - Tumor heterogeneity ####
library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
setwd("~/Dropbox/THCA_scrna/infercnv_epi_sub3/shannon/")
read.table("infercnv.observation_groupings_PTC_type1_k1.txt") -> PTC_type1_k1
read.table("infercnv.observation_groupings_PTC_type1_k2.txt") -> PTC_type1_k2
read.table("infercnv.observation_groupings_PTC_type1_k3.txt") -> PTC_type1_k3
read.table("infercnv.observation_groupings_PTC_type1_k4.txt") -> PTC_type1_k4
read.table("infercnv.observation_groupings_PTC_type2_k1.txt") -> PTC_type2_k1
read.table("infercnv.observation_groupings_PTC_type2_k2.txt") -> PTC_type2_k2
read.table("infercnv.observation_groupings_ATC_Tumor_k1.txt") -> ATC_k1
read.table("infercnv.observation_groupings_ATC_Tumor_k2.txt") -> ATC_k2
read.table("infercnv.observation_groupings_ATC_Tumor_k3.txt") -> ATC_k3
read.table("infercnv.observation_groupings_ATC_Tumor_k4.txt") -> ATC_k4

ATC_k1.2 <- rbind(ATC_k1,ATC_k2)

PTC_type1_k1[,2] <- "type1_k1"
PTC_type1_k2[,2] <- "type1_k2"
PTC_type1_k3[,2] <- "type1_k3"
PTC_type1_k4[,2] <- "type1_k4"
PTC_type2_k1[,2] <- "type2_k1"
PTC_type2_k2[,2] <- "type2_k2"
ATC_k1.2[,2] <- "ATC_k1.2"
ATC_k3[,2] <- "ATC_k3"
ATC_k4[,2] <- "ATC_k4"

total <- rbind(PTC_type1_k1, PTC_type1_k2, PTC_type1_k3, PTC_type1_k4, PTC_type2_k1, PTC_type2_k2, ATC_k1.2, ATC_k3, ATC_k4)
colnames(total) <- c("cell", "dendrogram")
row.names(total) <- total$cell
total[,1] <- NULL

Idents(Mye_thy) <- 'orig.ident'
THCA_only <- subset(Mye_thy, subset = infercnv_group %in% c("PTC type1", "PTC type2", "ATC"), idents = c('PTC2', 'PTC3', 'PTC4', 'ATC1', 'ATC2', 'ATC3'))
THCA_only <- AddMetaData(THCA_only, total)
THCA_only <- SetIdent(THCA_only, value="dendrogram")
THCA_only_subtypes_per_donor <- THCA_only@meta.data %>% select(orig.ident, dendrogram) %>% group_by(orig.ident, dendrogram) %>% summarise(n = n()) %>% as.data.frame()
THCA_only_subtypes_per_donor <- THCA_only_subtypes_per_donor %>% group_by(orig.ident) %>% mutate(prop = n / sum(n)) %>% as.data.frame()
THCA_only_subtypes_per_donor$prop <- round(THCA_only_subtypes_per_donor$prop, 3)
THCA_only_subtypes_per_donor <- plyr::ddply(THCA_only_subtypes_per_donor, "orig.ident", transform, label_y=cumsum(prop))
THCA_only_subtypes_per_donor$dendrogram <- factor(THCA_only_subtypes_per_donor$dendrogram, 
                                                  levels = c("type1_k1","type1_k2","type1_k3","type1_k4","type2_k1","type2_k2","ATC_k1", "ATC_k2", "ATC_k3", "ATC_k4"))
library(openxlsx)
write.xlsx(THCA_only_subtypes_per_donor, file = '/Users/yoo/Dropbox/THCA_scrna/infercnv_epi_sub3/shannon/20240428_THCA_PTC_ATC_subtypes_per_donor.xlsx')
THCA_only_subtypes_per_donor <- read.xlsx('/Users/yoo/Dropbox/THCA_scrna/infercnv_epi_sub3/shannon/20240428_THCA_PTC_ATC_subtypes_per_donor.xlsx')

library(vegan)
shannon <- c()
patient <- c()
for (i in unique(THCA_only_subtypes_per_donor$orig.ident)) {
  shannon <- c(shannon, round(diversity(THCA_only_subtypes_per_donor[THCA_only_subtypes_per_donor$orig.ident == i,]$n), 3))
  patient <- c(patient, i)
}
shannon_patient <- data.frame(patient, shannon)
shannon_order <- order(shannon, decreasing = T)
shannon_patient$patient[shannon_order]
shannon_patient$patient <- factor(shannon_patient$patient, levels = c(shannon_patient$patient[shannon_order]))

ggplot(shannon_patient) + 
  geom_bar(stat = "identity", aes(x=patient, y = shannon, fill=patient), color = "black", size = 0.3) +
  theme_classic() + scale_fill_manual(values = c("#F8766D", "#9590FF", "#C09B00", "#00BF7D", "#00BAE0", "#E76BF3")) +
  theme(axis.title.x = element_text(size=28, vjust = 0, family="Arial", margin = margin(0.28, 0, 0, 0, "cm")), 
        axis.text.x = element_text(family="Arial", size = 23, angle = 45, hjust = 1, vjust = 1, color = 'black'), 
        axis.title.y = element_text(size=28, vjust = 0, family="Arial", margin = margin(0, 0.28, 0, 0, "cm")), 
        axis.text.y = element_text(family="Arial", size = 23, hjust = 1, vjust = 1, color = 'black'), panel.grid = element_blank(),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", size=23),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.key.size = unit(0.4, 'cm')) + 
  labs(x = "Patient", y = "Shannon diversity index") + 
  geom_text(aes(x=patient, y=shannon, label=shannon), colour="black", family="Arial", size = 8, position = position_stack(vjust=0.5))+
  guides(fill = guide_legend(override.aes = list(colour = "black", size = 5), ncol = 6)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4), labels= c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4"))





