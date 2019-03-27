#### FIGURE 2 ####

#load libraries
pacman::p_load(rio, tidyverse, pheatmap, data.table, Seurat, readxl, GOfuncR)

#custom functions
#--------------------------Read data------------------------------
readat <- function(id) {
  pid = paste0("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Figure.2/", id)
  dat = import(pid)
  dat = dat %>% remove_rownames %>% column_to_rownames(var = names(dat)[1])
}

readex <- function(id, sheet) {
  pid = paste0("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Figure.2/", id)
  dat = read_excel(pid, sheet = sheet)
}

#----------------------------Theme1--------------------------------
theme_jo <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      legend.key = element_blank(), 
      strip.background = element_blank(),
      strip.text = element_text(size = 20, color = "grey10"),
      plot.title = element_blank(), 
      axis.title = element_text(size = 30, color = "grey10"), 
      axis.text = element_text(size = 25, color = "grey10"), 
      legend.title = element_blank(), 
      legend.text = element_blank(),
      axis.line = element_line(size = 1.1, color = "grey10"),
      axis.ticks = element_line(size = 1.1, color = "grey10"),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none")
}

#----------------------------Theme2--------------------------------
theme_blank <- function(base_size = 12){
    theme_jo(base_size = base_size) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank()
      )
}

#----------------------------Darken--------------------------------
darken <- function(color, factor=1.5){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) <- names(color)
  col
}

#### FIGURE 2A ####

#read count data
dat.ref.10x.count = Read10X("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Additions/TFi.Ci.10x.g19")

#seurat
iNS.10x = CreateSeuratObject(raw.data = dat.ref.10x.count, 
                             min.cells = 3, 
                             min.genes = 200, 
                             normalization.method = "LogNormalize", 
                             display.progress = T)

#regress on mitochondrial genes and number of UMI
mito.genes = grep(pattern = "MT-", x = rownames(x = iNS.10x@data), value = TRUE)
percent.mito = Matrix::colSums(iNS.10x@raw.data[mito.genes, ]) / Matrix::colSums(iNS.10x@raw.data)

#cluster cells
iNS.10x = AddMetaData(object = iNS.10x, metadata = percent.mito, col.name = "percent.mito") %>%
  FilterCells(subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(4500, 0.2)) %>%
  FindVariableGenes(x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 1, do.plot = F, display.progress = T) %>%
  ScaleData(genes.use = .@var.genes, vars.to.regress = c("nUMI", "percent.mito"), display.progress = T) %>%
  RunPCA(pcs.compute = 30, do.print = F, pcs.print = 0, genes.print = 0) %>%
  ProjectPCA(do.print = FALSE, pcs.print = 0, pcs.store = 30, genes.print = 0) %>%
  FindClusters(dims.use = 1:20, resolution = 0.7, print.output = T, force.recalc = T) %>%
  RunTSNE(dims.use = 1:20, do.fast = T)
#TSNEPlot(iNS.10x, do.label = T)

#reassigne cluster identities
iNS.10x@meta.data$res.0.7 = chartr("021368574", "123456789", iNS.10x@meta.data$res.0.7)
new.ident = as.factor(chartr("021368574", "123456789", iNS.10x@ident))
names(new.ident) = names(iNS.10x@ident)
iNS.10x@ident = new.ident

#create data frame for ggplot
iNS.10x.df = as.data.frame(iNS.10x@dr$tsne@cell.embeddings) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(unlist(lapply(strsplit(Cell, "[.]"), function(X) X[[1]])), levels = c("TFi", "Ci")),
         Cluster = factor(iNS.10x@meta.data$res.0.7, levels = c(1,2,3,4,5,6,7,8,9)))

#create color palettes
my.pal.10x.tsne.fill = c("1" = "#063c8d", "2" = "#0a8873", "3" = "#810081", "4" = "#c51385", "5" = "#00919e", 
                         "6" = "#00b900", "7" = "#f55701", "8" = "#f4ea0b", "9" = "#df0909") 
my.pal.10x.tsne.col = darken(my.pal.10x.tsne.fill)

#plot tSNE
ggplot(iNS.10x.df, aes(tSNE_1, tSNE_2, fill = Cluster, color = Cluster, 
                       shape = as.factor(Sample))) +
  geom_point(alpha = 1, size = 3.5) +
  scale_fill_manual(values = my.pal.10x.tsne.fill) +
  scale_color_manual(values = my.pal.10x.tsne.col) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(5, 4)) +
  theme_jo()

#### FIGURE 2B ####

#read data
dat.ref.10x = readat("dat.reference.10x.qc.cpm.log2.csv")

#define genes
marker.g.2 = c("NRCAM", "SYT1", "SERPINI1", "SST", "SYNGR1", "SYNPO2", "SPOCK1", "STMN2",
               "SNAI2", "COL15A1", "CPXM1", "COL1A1", "EDIL3", "HAPLN1", "THY1", "MFAP4")

#create matrix
dat.marker.2 = dat.ref.10x[marker.g.2, ] %>%
  t %>%
  cbind(iNS.10x.df, .) %>%
  arrange(tSNE_2) %>%
  mutate(cut = cut(.$tSNE_2, 40)) %>%
  split(as.factor(.$cut)) %>%
  lapply(., function(X) colMeans(X[, marker.g.2])) %>%
  do.call("rbind", .)
dat.marker.2 = dat.marker.2 %>%
  .[rev(row.names(.)), ]

#plot heatmap
pheatmap(dat.marker.2, 
         cluster_row = F, 
         scale = "column", 
         border_color = "grey60", 
         show_rownames = F, 
         fontsize = 13,
         cutree_cols = 2,
         color = colorRampPalette(c("#001532", "#001532","white", "#e99400", "#e99400"))(100), 
         legend = T)

#### FIGURE 2C ####

#identify cluster-specific markers
cluster.markers = FindAllMarkers(iNS.10x, only.pos = T)
cluster.markers$p_val_log = -log10(cluster.markers$p_val)

#create matrix for heatmap
top10 = cluster.markers %>% group_by(cluster) %>% top_n(30, p_val_log)
dat.heat = DoHeatmap(object = iNS.10x,
                     genes.use = top10$gene,
                     slim.col.label = TRUE,
                     remove.key = F)

dat.heat.df = dat.heat$data[, -4]
dat.heat.df = dcast(dat.heat.df, gene~cell) %>%
  column_to_rownames(var = "gene")
gene.ident = top10 %>%
  dplyr::select(c(gene, cluster, p_val_log))
cell.ident = iNS.10x.df %>%
  dplyr::select(c(Cell, Cluster, Sample)) %>%
  mutate(Ident = as.numeric(as.character(Cluster))) %>%
  arrange(Ident)
m = data.frame(gene = row.names(dat.heat.df)) %>%
  left_join(gene.ident, by = "gene") %>%
  arrange(cluster, desc(p_val_log))
dat.heat.df = dat.heat.df[m$gene, cell.ident$Cell]
pheatmap(as.matrix(dat.heat.df), cluster_cols = F, cluster_rows = F)
dat.ann = data.frame(Sample = cell.ident$Sample,
                     Cluster = chartr("123456789", "ABCDEFGHI", cell.ident$Cluster))
row.names(dat.ann) = cell.ident$Cell

#create color annotations
ann_colors = list(
  Cluster = c("A" = "#063c8d", "B" = "#0a8873", "C" = "#810081", "D" = "#c51385", "E" = "#00919e",
              "F" = "#00b900", "G" = "#f55701", "H" = "#f4ea0b", "I" = "#df0909"),
  Sample = c("Ci" = "#003829", "TFi" = "#ccc200"))

#plot heatmap
pheatmap(as.matrix(dat.heat.df), 
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_col = dat.ann, 
         annotation_colors = ann_colors,
         annotation_legend = F,
         color = colorRampPalette(c("white", "#728fbd", "#001532", "#ffd58f", "#e99400"))(100),
         show_rownames = T,
         show_colnames = F,
         fontsize = 15,
         fontsize_row = 4,
         gaps_col = c(1316, 2014, 2994, 3331, 3443, 3471, 3592, 3646),
         gaps_row = c(4, 16, 35, 59, 80, 98, 125, 145))