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
########################################## FIGURE S2A-C ############################################
####################################################################################################

#read 10x data, create Seurat
dat_ref_10x_count <- Read10X(file.path(my_path, "Data", "TFi_Ci_10x_g19"))

#seurat
iNS_10x <- CreateSeuratObject(raw.data = dat_ref_10x_count, 
                              min.cells = 3, 
                              min.genes = 200, 
                              normalization.method = "LogNormalize", 
                              display.progress = F)

#regress on mitochondrial genes and number of UMI
mito_genes <- grep(pattern = "MT-", x = rownames(x = iNS_10x@data), value = TRUE)
percent_mito <- Matrix::colSums(iNS_10x@raw.data[mito_genes, ]) / Matrix::colSums(iNS_10x@raw.data)
iNS_10x <- AddMetaData(object = iNS_10x, metadata = percent_mito, col.name = "percent.mito")

#get 10x qc data
dat_viol <- VlnPlot(iNS_10x, c("nGene", "nUMI", "percent.mito"), nCol = 3, do.return = T, return.plotlist = T)

#extract gene data
dat_gene <- dat_viol[[1]]$data %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(c(rep("Ci", times = length(grep("Ci", .$Cell))), 
                           rep("TFi", times = length(grep("TFi", .$Cell)))), levels = c("Ci", "TFi")))

#extract umi data
dat_umi <- dat_viol[[2]]$data %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(c(rep("Ci", times = length(grep("Ci", .$Cell))), 
                           rep("TFi", times = length(grep("TFi", .$Cell)))), levels = c("Ci", "TFi")))

#extract mito data
dat_mito <- dat_viol[[3]]$data %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(c(rep("Ci", times = length(grep("Ci", .$Cell))), 
                           rep("TFi", times = length(grep("TFi", .$Cell)))), levels = c("Ci", "TFi")))

#color palette
qc_pal_fill <- c("#6d6d6d", "#c0282c")
qc_pal_col <- c("#202020", "#73171b")

#plot gene data
ggplot(dat_gene, aes(Sample, feature, color = Sample, fill = Sample)) +
  scale_fill_manual(values = qc_pal_fill) +
  scale_color_manual(values = qc_pal_col) +
  geom_jitter(shape = 21, size = 2, alpha = 1, width = 0.3) +
  geom_violin(scale = "count", lwd = 0.8, fill = "grey90", alpha = 0.1) +
  geom_hline(yintercept = 4500, color = "grey10", size = 1, lty = 2) +
  geom_hline(yintercept = 1000, color="grey10", size = 1, lty = 2) +
  theme_jo()

#plot umi data
ggplot(dat_umi, aes(Sample, feature, color = Sample, fill = Sample)) +
  scale_fill_manual(values = qc_pal_fill) +
  scale_color_manual(values = qc_pal_col) +
  geom_jitter(shape = 21, size = 2, alpha = 1, width = 0.3) +
  geom_violin(scale = "count", lwd = 0.8, fill = "grey90", alpha = 0.1) +
  theme_jo()

#plot mito data
ggplot(dat_mito, aes(Sample, feature, color = Sample, fill = Sample)) +
  scale_fill_manual(values = qc_pal_fill) +
  scale_color_manual(values = qc_pal_col) +
  geom_jitter(shape = 21, size = 2, alpha = 1, width = 0.3) +
  geom_violin(scale = "count", lwd = 0.8, fill = "grey90", alpha = 0.1) +
  geom_hline(yintercept = 0.2, color = "grey10", size = 1, lty = 2)+
  theme_jo()

####################################################################################################
########################################## FIGURE S2A-C ############################################
####################################################################################################

#apply thresholds, find variable genes
iNS_10x <- iNS_10x %>%
  FilterCells(subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(4500, 0.2)) %>%
  NormalizeData(display.progress = F) %>%
  FindVariableGenes(x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 1, do.plot = F, display.progress = F) 

#create data frame for plotting
dat_var <- iNS_10x@hvg.info %>%
  mutate(label = row.names(.))

#plot
ggplot(dat_var, aes(gene.mean, gene.dispersion.scaled, label = label)) +
  geom_point(size = 3, fill = "grey70", color = "grey30", shape = 21, alpha = 1) +
  geom_point(data = dat_var[which(dat_var$gene.mean >= 0.01 & dat_var$gene.mean <= 8 & dat_var$gene.dispersion.scaled >= 1), ], 
             size = 3, fill = "#eb9926", color = "#a1691a", shape = 21, alpha = 1) +
  geom_text(data = dat_var[which(dat_var$gene.mean > 1.5 & dat_var$gene.dispersion.scaled > 3), ], size = 5, color = "grey20", 
            alpha = 1, nudge_x = 0.3, nudge_y = 0.5) +
  geom_text(data = dat_var[which(dat_var$gene.mean < 1.6 & dat_var$gene.mean > 0.3 & dat_var$gene.dispersion.scaled > 6), ], size = 5, color = "grey10", 
            alpha = 1, nudge_x = 0.3, nudge_y = 0.5) +
  geom_hline(yintercept = 1, color="grey10", size = 1, lty = 2) +
  geom_vline(xintercept = 0.01, color = "grey10", size = 1, lty = 2) +
  theme_jo()