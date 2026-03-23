setwd('/Users/yoo/Dropbox/THCA_scrna/')
library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(scales)

#Fig. 1B - global UMAP####

THCA <- readRDS('thyroid_scrnaseq/processed.rds')

DimPlot(THCA, group.by = "annotation_v5") & 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12.5, family="Arial"),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) &
  guides(col=guide_legend(ncol=3, override.aes = list(size = 3))) &
  NoAxes() &
  ggtitle(NULL)


# Fig. 1C - cell type marker genes ####
THCA <- readRDS('thyroid_scrnaseq/processed.rds')

folli <- c('TG', 'TPO', 'TFF3')
PTC <- c('NPC2','GDF15','SLC34A2')
ATC <- c('UBE2C','S100A10','TFPI2')
myeloid <- c('LYZ', 'CD14', 'S100A9')
T.NK <- c('CD3D','CD3E','CD3G','NKG7')
fibro <- c('COL1A1','COL1A2','ACTA2')
B_cell <- c('MS4A1','BANK1','VPREB3','CD79A')
plasma <- c('IGLL5','MZB1','IGJ')
prolif <- c('HMGB2','HMGN2','STMN1')
endo <- c('VWF','PLVAP','CLDN5')
pDC <- c('LILRA4','GZMB','IL3RA')

markers <- c(folli,PTC,ATC,T.NK,B_cell,plasma,myeloid,pDC,fibro,endo,prolif)

Idents(THCA) <- 'annotation_v5'
DotPlot(THCA, features = markers, group.by = "annotation_v5") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 12)) +
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "Reds"))


#Fig. 1D - cell type freq. changes across samples  ####
THCA <- readRDS('thyroid_scrnaseq/processed.rds')

Idents(THCA) <- 'orig.ident'
table_THCA <- THCA@meta.data %>% dplyr::select(orig.ident, annotation_v5) %>% table()
table_THCA_v2 <- ( table_THCA / table_THCA %>% apply(1,sum) ) %>% data.frame()
table_THCA_v2$orig.ident <- factor(table_THCA_v2$orig.ident, levels = c("PTC3", "PTC1", "PTC4", "PTC2", "PTCEO", "ATC2", "ATC1", "ATC3"))

#barplot
ggplot(table_THCA_v2, aes(fill=annotation_v5, y=Freq, x=orig.ident)) + 
  geom_bar(position="fill", stat="identity") +
  theme(plot.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 30, size = 15, vjust = 0.9, hjust = 0.95),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


#Fig. 1E - cell type freq. changes across cancer types ####
THCA <- readRDS('thyroid_scrnaseq/processed.rds')

Idents(THCA) <- 'orig.ident'
table_THCA <- THCA@meta.data %>% select(orig.ident, annotation_v5) %>% table()
table_THCA_v2 <- ( table_THCA / table_THCA %>% apply(1,sum) ) %>% data.frame()
table_THCA_v2 <- table_THCA_v2 %>% mutate(group = case_when(
  grepl('^PTC',table_THCA_v2$orig.ident) == TRUE ~ 'PTC',
  grepl('^ATC',table_THCA_v2$orig.ident) == TRUE ~ 'ATC'
))
table_THCA_v2$group <- factor(table_THCA_v2$group, levels = c('PTC','ATC'))

ggplot(table_THCA_v2, aes(fill=annotation_v5, y=Freq, x=group)) + 
  geom_bar(position="fill", stat="identity") +
  guides(fill = guide_legend(ncol = 3)) +
  theme(plot.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 30, size = 15, vjust = 0.9, hjust = 0.95),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


#Fig. 1F - RCTD results####
ST_merge <- readRDS("R_objects/thca_visium_merge.rds")

weights_top_all <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("celltype"))
for (smp in unique(ST_merge$sample)) {
  myRCTD <- readRDS(paste0("visium_public/rctd/RCTD_full_", smp, ".rds"))
  weights <- myRCTD@results$weights %>% as.data.frame()
  weights_top <- data.frame(
    celltype = colnames(weights)[apply(weights, 1, which.max)],
    row.names = rownames(weights))
  weights_top_all <- rbind(weights_top_all,weights_top)
}

ST_merge <- AddMetaData(ST_merge, weights_top_all)

mycol <- scPalette(length(unique(ST_merge$celltype)))
names(mycol) <- unique(ST_merge$celltype)
show_col(mycol)

mycol['myCAF'] <- "#4385BF"
mycol['Plasma'] <- "#6B5F88"
mycol[7] <- '#999999'
mycol['IL-8+ Monocyte'] <- "#87638F"
mycol['Follicular cell'] <- "#53A651"
mycol['ATC'] <- "#E41A1C"
mycol['Mature naive B'] <- "#683A68"
mycol['PTC type1'] <- '#F29022'
mycol['PTC type2'] <- "#C26F57"
mycol['Classical Monocytes'] <- "#D690C6" 
mycol['SPP1+ TAM'] <- "#DDE029"
mycol['iCAF'] <- "#3F918D"
mycol['CD8+ proliferative T'] <- "#7A3F91"

SpatialDimPlot(ST_merge, group.by = 'celltype', images = names(ST_merge@images), ncol = 8, cols = mycol) & NoLegend()


#Fig. 1G - RCTD cell type freq. changes across cancer states####
library(rstatix)
library(ggplot2)
library(ggsignif)

ST_merge <- readRDS("R_objects/thca_visium_merge.rds")
Idents(ST_merge) <- 'sample'
table_ST_merge <- ST_merge@meta.data %>% select(sample, celltype) %>% table()
table_ST_merge_v2 <- ( table_ST_merge / table_ST_merge %>% apply(1,sum) ) %>% data.frame()
celltypes_flt <- table_ST_merge_v2 %>% filter(Freq > 0.1) %>% pull(celltype) %>% unique()
table_ST_merge_v2 <- table_ST_merge_v2 %>% filter(celltype %in% c('iCAF','myCAF'))
table_ST_merge_v2 <- table_ST_merge_v2 %>% mutate(group = case_when(
  grepl('^N',table_ST_merge_v2$sample) == TRUE ~ 'PT',
  grepl('^PTC',table_ST_merge_v2$sample) == TRUE ~ 'PTC',
  grepl('^LPTC',table_ST_merge_v2$sample) == TRUE ~ 'LPTC',
  grepl('^ATC',table_ST_merge_v2$sample) == TRUE ~ 'ATC'
))
table_ST_merge_v2$group <- factor(table_ST_merge_v2$group, levels = c('PT','PTC','LPTC','ATC'))

#Wilcoxon test
wilcox_compare <- function(df, celltype_name, group1, group2) {
  wilcox.test(
    Freq ~ group,
    data = df %>% 
      filter(celltype == celltype_name, group %in% c(group1, group2)),
    alternative = "less"  # PT < cancer group (Freq higher in cancer)
  )
}

table_ST_merge_v2 %>%
  group_by(celltype) %>%
  wilcox_test(
    Freq ~ group,
    ref.group = "PT",
    alternative = "less"
  ) %>%
  adjust_pvalue(method = "BH")

ggplot(data=table_ST_merge_v2,aes(x=factor(celltype), y=Freq, fill=group)) +
  geom_boxplot() +
  theme_classic() +
  labs(y = 'Fraction', fill = 'Cancer state') +
  theme(axis.text.x = element_text(colour = "black", angle = 45, size = 27, vjust=1, hjust=1), 
        axis.text.y = element_text(colour = "black", size = 27), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(colour = "black", size = 27),
        plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 27),
        legend.title = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        legend.spacing = unit(1, 'cm'), legend.spacing.x = unit(0.5,'cm')) +
  geom_signif(
    y_position = c(0.33,0.36), xmin = c(1.7,1.7), xmax = c(2.125, 1.875),
    annotation = c('P = 0.087', 'P = 0.087'), 
    tip_length = 0, textsize = 8, color = 'black', vjust=-0.1) +
  scale_fill_manual(values = c("PT" = "#4DAF4A",
                               "PTC" = "#377EB8",
                               "LPTC" = "#984EA3",
                               "ATC" = "#E41A1C")) +
  guides(
    fill  = guide_legend(ncol = 2)
  )

