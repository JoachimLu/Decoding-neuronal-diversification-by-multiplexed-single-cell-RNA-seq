#### FIGURE S2 ####

#load libraries
pacman::p_load(rio, tidyverse, data.table, Seurat, readxl)

#custom functions
#--------------------------Read data------------------------------
readat <- function(id) {
  pid = paste0("C:/Users/joachim/Desktop/scRNAseq.iN/Data/Figure.S2/", id)
  dat = import(pid)
  dat = dat %>% remove_rownames %>% column_to_rownames(var = names(dat)[1])
}

readex <- function(id, sheet) {
  pid = paste0("C:/Users/joachim/Desktop/scRNAseq.iN/Data/Figure.S2/", id)
  dat = read_excel(pid, sheet = sheet)
}

#----------------------------Theme1--------------------------------
theme_jo <- function(base_size = 12){
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
theme_blank = function(base_size = 12){
  theme_jo(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      strip.text = element_blank()
    )
}

#### FIGURE S2A-C ####

#read 10x data, create Seurat
dat.ref.10x.count = Read10X("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Additions/TFi.Ci.10x.g19")

#seurat
iNS.10x = CreateSeuratObject(raw.data = dat.ref.10x.count, 
                             min.cells = 3, 
                             min.genes = 200, 
                             normalization.method = "LogNormalize", 
                             display.progress = F)

#regress on mitochondrial genes and number of UMI
mito.genes = grep(pattern = "MT-", x = rownames(x = iNS.10x@data), value = TRUE)
percent.mito = Matrix::colSums(iNS.10x@raw.data[mito.genes, ]) / Matrix::colSums(iNS.10x@raw.data)
iNS.10x = AddMetaData(object = iNS.10x, metadata = percent.mito, col.name = "percent.mito")

#get 10x qc data
dat.viol = VlnPlot(iNS.10x, c("nGene", "nUMI", "percent.mito"), nCol = 3, do.return = T, return.plotlist = T)

#extract gene data
dat.gene = dat.viol[[1]]$data %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(c(rep("Ci", times = length(grep("Ci", .$Cell))), 
                           rep("TFi", times = length(grep("TFi", .$Cell)))), levels = c("Ci", "TFi")))

#extract umi data
dat.umi = dat.viol[[2]]$data %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(c(rep("Ci", times = length(grep("Ci", .$Cell))), 
                           rep("TFi", times = length(grep("TFi", .$Cell)))), levels = c("Ci", "TFi")))

#extract mito data
dat.mito = dat.viol[[3]]$data %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(c(rep("Ci", times = length(grep("Ci", .$Cell))), 
                           rep("TFi", times = length(grep("TFi", .$Cell)))), levels = c("Ci", "TFi")))

#color palette
qc.pal.fill = c("#6d6d6d", "#c0282c")
qc.pal.col = c("#202020", "#73171b")

#plot gene data
ggplot(dat.gene, aes(Sample, feature, color = Sample, fill = Sample)) +
  scale_fill_manual(values = qc.pal.fill) +
  scale_color_manual(values = qc.pal.col) +
  geom_jitter(shape = 21, size = 2, alpha = 1, width = 0.3) +
  geom_violin(scale = "count", lwd = 0.8, fill = "grey90", alpha = 0.1) +
  geom_hline(yintercept = 4500, color = "grey10", size = 1, lty = 2) +
  geom_hline(yintercept = 1000, color="grey10", size = 1, lty = 2) +
  theme_jo()

#plot umi data
ggplot(dat.umi, aes(Sample, feature, color = Sample, fill = Sample)) +
  scale_fill_manual(values = qc.pal.fill) +
  scale_color_manual(values = qc.pal.col) +
  geom_jitter(shape = 21, size = 2, alpha = 1, width = 0.3) +
  geom_violin(scale = "count", lwd = 0.8, fill = "grey90", alpha = 0.1) +
  theme_jo()

#plot mito data
ggplot(dat.mito, aes(Sample, feature, color = Sample, fill = Sample)) +
  scale_fill_manual(values = qc.pal.fill) +
  scale_color_manual(values = qc.pal.col) +
  geom_jitter(shape = 21, size = 2, alpha = 1, width = 0.3) +
  geom_violin(scale = "count", lwd = 0.8, fill = "grey90", alpha = 0.1) +
  geom_hline(yintercept = 0.2, color = "grey10", size = 1, lty = 2)+
  theme_jo()

#### FIGURE S2D ####

#apply thresholds, find variable genes
iNS.10x = iNS.10x %>%
  FilterCells(subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(4500, 0.2)) %>%
  NormalizeData(display.progress = F) %>%
  FindVariableGenes(x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 1, do.plot = F, display.progress = F) 

#create data frame for plotting
dat.var = iNS.10x@hvg.info %>%
  mutate(label = row.names(.))

#plot
ggplot(dat.var, aes(gene.mean, gene.dispersion.scaled, label = label)) +
  geom_point(size = 3, fill = "grey70", color = "grey30", shape = 21, alpha = 1) +
  geom_point(data = dat.var[which(dat.var$gene.mean >= 0.01 & dat.var$gene.mean <= 8 & dat.var$gene.dispersion.scaled >= 1), ], 
             size = 3, fill = "#eb9926", color = "#a1691a", shape = 21, alpha = 1) +
  geom_text(data = dat.var[which(dat.var$gene.mean > 1.5 & dat.var$gene.dispersion.scaled > 3), ], size = 5, color = "grey20", 
            alpha = 1, nudge_x = 0.3, nudge_y = 0.5) +
  geom_text(data = dat.var[which(dat.var$gene.mean < 1.6 & dat.var$gene.mean > 0.3 & dat.var$gene.dispersion.scaled > 6), ], size = 5, color = "grey10", 
            alpha = 1, nudge_x = 0.3, nudge_y = 0.5) +
  geom_hline(yintercept = 1, color="grey10", size = 1, lty = 2) +
  geom_vline(xintercept = 0.01, color = "grey10", size = 1, lty = 2) +
  theme_jo()