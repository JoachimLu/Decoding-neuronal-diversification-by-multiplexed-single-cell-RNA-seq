#check if pacman package is loaded and install if not
if (!require(pacman)) {
  install.packages("pacman")
  if(!require(pacman)) stop("Package pacman not found")
}

#set data path
my_path <- "D:/Projects/Single_Cell_Convert_Seq"

#load libraries and custom functions
pacman::p_load(rio, tidyverse, pheatmap, data.table, Seurat, readxl, GOfuncR, Biobase)
source(file.path(my_path, "R", "utilities.R"))

#check package version of Seurat
if(package.version("Seurat") != "2.3.4") {
  stop("Analysis in this script is based on Seurat version 2.3.4 and wont run on newer versions")
}

####################################################################################################
############################################ FIGURE 2A #############################################
####################################################################################################

#read Seurat object
iNS_10x <- get(load(file.path(my_path, "Data", "Seurat_iN_10x.Rmd")))

#create data frame for ggplot
iNS_10x_df <- as.data.frame(iNS_10x@dr$umap@cell.embeddings) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(unlist(lapply(strsplit(Cell, "[.]"), function(X) X[[1]])), levels = c("Ci", "TFi")),
         Cluster = factor(iNS_10x@meta.data$res.0.7, levels = c(1,2,3,4,5,6,7,8,9)))

#create color palettes
my_pal_10x_tsne_fill <- c("1" = "#5e5043", 
                          "2" = "#d6c8a5", 
                          "3" = "#810081", 
                          "4" = "#c51385", 
                          "5" = "#00b900", 
                          "6" = "#f55701", 
                          "7" = "#27588a", 
                          "8" = "#f4c00a", 
                          "9" = "#df0909") 
my_pal_10x_tsne_col <- darken(my_pal_10x_tsne_fill)

#plot tSNE
ggplot(iNS_10x_df, aes(UMAP1, UMAP2, fill = Cluster, color = Cluster, 
                       shape = as.factor(Sample))) +
  geom_point(alpha = 1, size = 3.5) +
  scale_fill_manual(values = my_pal_10x_tsne_fill) +
  scale_color_manual(values = my_pal_10x_tsne_col) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(5, 4)) +
  theme_jo()

####################################################################################################
############################################ FIGURE 2B #############################################
####################################################################################################

#read data
dat_ref_10x <- readat(file.path(my_path, "Data", "dat_reference_10x_qc_cpm_log2.csv"))

#define genes
marker_g <- c("NRCAM", "SYT1", "SERPINI1", "SST", "SYNGR1", "SYNPO2", "SPOCK1",
              "SNAI2", "COL15A1", "CPXM1", "FST", "COL1A1", "CDH11", "VCAN")

#create matrix
dat_marker <- dat_ref_10x[marker_g, ] %>%
  t %>%
  cbind(iNS_10x_df, .) %>%
  arrange(UMAP1) %>%
  mutate(cut = cut(.$UMAP1, 40)) %>%
  split(as.factor(.$cut)) %>%
  lapply(., function(X) colMeans(X[, marker_g])) %>%
  do.call("rbind", .) %>%
  t

#plot heatmap
pheatmap(dat_marker, 
         cluster_col = F, 
         scale = "row", 
         border_color = "grey60", 
         show_colnames = F, 
         fontsize = 13,
         cutree_rows = 2,
         color = colorRampPalette(c("#cdcdcb", "#cdcdcb", "#e7f5b6", "#3db2c4", "#0d3286", "#0d3286"))(100), 
         legend = T)

####################################################################################################
############################################ FIGURE 2C #############################################
####################################################################################################

#identify cluster-specific markers
cluster_markers <- FindAllMarkers(iNS_10x, only.pos = T)
cluster_markers$p_val_log = -log10(cluster_markers$p_val)

#create matrix for heatmap
top10 <- cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(30, p_val_log)
dat_heat <- DoHeatmap(object = iNS_10x,
                      genes.use = top10$gene,
                      slim.col.label = TRUE,
                      remove.key = F)

dat_heat_df <- dat_heat$data[, -4]
dat_heat_df <- dcast(dat_heat_df, gene~cell) %>%
  column_to_rownames(var = "gene")
gene_ident <- top10 %>%
  dplyr::select(c(gene, cluster, p_val_log))
cell_ident <- iNS_10x_df %>%
  dplyr::select(c(Cell, Cluster, Sample)) %>%
  mutate(Ident = as.numeric(as.character(Cluster))) %>%
  arrange(Ident)
m <- data.frame(gene = row.names(dat_heat_df)) %>%
  left_join(gene_ident, by = "gene") %>%
  arrange(cluster, desc(p_val_log))
dat_heat_df <- dat_heat_df[m$gene, cell_ident$Cell]
pheatmap(as.matrix(dat_heat_df), cluster_cols = F, cluster_rows = F)
dat_ann <- data.frame(Sample = cell_ident$Sample,
                      Cluster = chartr("123456789", "ABCDEFGHI", cell_ident$Cluster))
row.names(dat_ann) <- cell_ident$Cell

#create color annotations
ann_colors <- list(
  Cluster = c("A" = "#5e5043", "B" = "#d6c8a5", "C" = "#810081", "D" = "#c51385", "E" = "#00919e",
              "F" = "#00b900", "G" = "#f55701", "H" = "#f4ea0b", "I" = "#df0909"),
  Sample = c("Ci" = "#415934", "TFi" = "#bacd79"))

#plot heatmap
pheatmap(as.matrix(dat_heat_df), 
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_col = dat_ann, 
         annotation_colors = ann_colors,
         annotation_legend = F,
         color = colorRampPalette(c("#cdcdcb", "#e7f5b6", "#3db2c4", "#0d3286"))(100),
         show_rownames = T,
         show_colnames = F,
         fontsize = 15,
         fontsize_row = 4,
         gaps_col = c(1316, 2014, 2994, 3331, 3443, 3471, 3592, 3646),
         gaps_row = c(4, 16, 35, 59, 80, 98, 125, 145))

####################################################################################################
############################################ FIGURE 2D #############################################
####################################################################################################

#define marker gene
marker_genes_10x <- c("SPOCK1", "FRZB", "MKI67", "CENPA", "DRD4", "SEMA3G",
                      "CHRNA5", "CP", "CHD7", "DLL3", "ANK3", "NKAIN4")

dat_all <- matrix(nrow = 0, ncol = 7)
for(i in 1:length(marker_genes_10x)) {
  #create data frame
  dat_marker = iNS_10x_df %>%
    mutate(Gene = marker_genes_10x[i],
           TPM = as.numeric(dat_ref_10x[marker_genes_10x[i], ])) 
  dat_all <- rbind(dat_all, dat_marker)
}
dat_all$Gene <- factor(dat_all$Gene, levels = marker_genes_10x)

#plot
ggplot(dat_all, aes(Cluster, TPM, fill = Cluster, color = Cluster)) +
  geom_violin(scale = "width", lwd = 0.8, alpha = 1, trim = T) +
  scale_color_manual(values = my_pal_10x_tsne_col) +
  scale_fill_manual(values = my_pal_10x_tsne_fill) +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(dat_marker$Cluster))) +
  theme_jo() +
  theme(axis.line = element_line(size = 0.5), 
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_blank(),
        axis.title = element_blank()) +
  facet_wrap(~Gene, nrow = 1)

####################################################################################################
############################################ FIGURE 2E #############################################
####################################################################################################

#plot
ggplot() +
  geom_point(data = iNS_10x_df[which(iNS_10x_df$Cluster %in% c(1,2,3)), ], aes(UMAP1, UMAP2,
                                                                               shape = as.factor(Sample)),
             fill = "white",
             color = "grey50",
             size = 3.5) +
  geom_point(data = iNS_10x_df[which(iNS_10x_df$Cluster %in% c(4,5,6,7,8,9)), ], aes(UMAP1, UMAP2,
                                                                                     fill = as.factor(Cluster),
                                                                                     color = as.factor(Cluster),
                                                                                     shape = as.factor(Sample)),
             size = 3.5) +
  scale_fill_manual(values = my_pal_10x_tsne_fill) +
  scale_color_manual(values = my_pal_10x_tsne_col) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(5, 4)) +
  theme_jo()