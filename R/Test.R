#libraries
pacman::p_load(rio, tidyverse, pheatmap, data.table, RColorBrewer, Seurat, readxl, monocle, GOfuncR)

#read 10x data, create Seurat
dat.ref.10x.count <- Read10X("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Additions/TFi.Ci.10x.g19")

#seurat
iNS.10x <- CreateSeuratObject(raw.data = dat.ref.10x.count,
                              min.cells = 3,
                              min.genes = 200,
                              normalization.method = "LogNormalize",
                              display.progress = T)

#regress on mitochondrial genes and number of UMI
mito.genes <- grep(pattern = "MT-", x = rownames(x = iNS.10x@data), value = TRUE)
percent.mito <- Matrix::colSums(iNS.10x@raw.data[mito.genes, ]) / Matrix::colSums(iNS.10x@raw.data)

#cluster cells
iNS.10x <- AddMetaData(object = iNS.10x, metadata = percent.mito, col.name = "percent.mito") %>%
  FilterCells(subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(4500, 0.2)) %>%
  FindVariableGenes(x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 1, do.plot = F, display.progress = F) %>%
  ScaleData(genes.use = .@var.genes, vars.to.regress = c("nUMI", "percent.mito"), display.progress = F) %>%
  RunPCA(pcs.compute = 30, do.print = F, pcs.print = 0, genes.print = 0) %>%
  ProjectPCA(do.print = FALSE, pcs.print = 0, pcs.store = 30, genes.print = 0) %>%
  FindClusters(dims.use = 1:20, resolution = 0.7, print.output = F, force.recalc = T) %>%
  RunTSNE(dims.use = 1:20, do.fast = T)

#reassigne cluster identities
iNS.10x@meta.data$res.0.7 <- chartr("021368574", "123456789", iNS.10x@meta.data$res.0.7)
iNS.10x@ident <- setNames(as.factor(chartr("021368574", "123456789", iNS.10x@ident)), names(iNS.10x@ident))

#create data frame for ggplot
iNS.10x.df <- as.data.frame(iNS.10x@dr$tsne@cell.embeddings) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(unlist(lapply(strsplit(Cell, "[.]"), function(X) X[[1]])), levels = c("TFi", "Ci")),
         Cluster = factor(iNS.10x@meta.data$res.0.7, levels = c(1,2,3,4,5,6,7,8,9)))

#create color palettes
my.pal.10x.tsne.fill <- c("1" = "#063c8d", "2" = "#0a8873", "3" = "#810081", "4" = "#c51385", "5" = "#00919e", 
                         "6" = "#00b900", "7" = "#f55701", "8" = "#f4ea0b", "9" = "#df0909") 
my.pal.10x.tsne.col <- darken(my.pal.10x.tsne.fill)

#plot tSNE
ggplot(iNS.10x.df, aes(tSNE_1, tSNE_2, fill = Cluster, color = Cluster, 
                       shape = as.factor(Sample))) +
  geom_point(alpha = 1, size = 3.5) +
  scale_fill_manual(values = my.pal.10x.tsne.fill) +
  scale_color_manual(values = my.pal.10x.tsne.col) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(5, 4)) +
  theme_jo()