#### FIGURE S4 ####

#load libraries
pacman::p_load(rio, tidyverse, data.table, Seurat, readxl, dendextend)

#custom functions
#--------------------------Read data------------------------------
readat <- function(id) {
  pid = paste0("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Figure.S4/", id)
  dat = import(pid)
  dat = dat %>% remove_rownames %>% column_to_rownames(var = names(dat)[1])
}

readex <- function(id, sheet) {
  pid = paste0("C:/Users/joachim/Desktop/scRNAseq.iN/Data/Figure.S4/", id)
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

#----------------------------Multiple replacements--------------------------------

gsub2 <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, fixed = T)
  }
  result
}

#### FIGURE S4A ####

#read tc data
my.files = c("tc.mat.count.raw.csv", "tc.mat.tpm.raw.csv", "ercc.counts.raw.csv", "ercc.tpm.raw.csv", "tc.qc.low.csv")
dat.all = pbapply::pblapply(my.files, readat)
names(dat.all) = gsub("\\.csv$", "", my.files)

#remove dead cells and cell doublets from count data
dat.count = dat.all[["tc.mat.count.raw"]][ ,-(which(colnames(dat.all[["tc.mat.count.raw"]]) %in% dat.all[["tc.qc.low"]]$LowQC))]

#remove cells with low number of detected genes
gene.count = apply(dat.count, 2, function(y) length(which(y > 0)))
dat.all[["tc.counts.qc"]] = as.data.frame(dat.count[, -(which(gene.count < (mean(gene.count) - 3*sd(gene.count))))])
cells.qc = names(dat.all[["tc.counts.qc"]])
dat.all[["tc.tpm.qc"]] = dat.all[["tc.mat.tpm.raw"]][, which(names(dat.all[["tc.mat.tpm.raw"]]) %in% cells.qc)]
dat.all[["ercc.counts.qc"]] = dat.all[["ercc.counts.raw"]][, which(names(dat.all[["ercc.counts.raw"]]) %in% cells.qc)]
dat.all[["ercc.tpm.qc"]] = dat.all[["ercc.tpm.raw"]][, which(names(dat.all[["ercc.tpm.raw"]]) %in% cells.qc)]
dat.all[["tc.tpm.qc.log2"]] = log2(dat.all[["tc.tpm.qc"]] + 1)
dat.all[["ercc.tpm.qc.log2"]] = log2(dat.all[["ercc.tpm.qc"]] + 1)

#calculate percentage of ERCC spikes
dat.ratio = as.data.frame((100*(colSums(dat.all[["ercc.counts.qc"]])) / 
                             (colSums(dat.all[["ercc.counts.qc"]]) + colSums(dat.all[["tc.counts.qc"]])))) %>%
  mutate(Group = unlist(lapply(strsplit(row.names(.), "[.]"), function(X) X[[1]]))) %>%
  rownames_to_column(var="Cell")
names(dat.ratio) = c("Cell", "Ratio", "Group")

#melt data frame
dat.ratio.m = melt(dat.ratio, id = c("Ratio", "Group")) %>%
  mutate(Group = factor(.$Group, levels = c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4")))

#create color palette
my.pal.qc = c("FIB" = "#033782", "Ci1" = "#8b008b", "Ci2" = "#d40b91", "TFi1" = "#ec1717", 
              "TFi2" = "#fa7600", "TFi3" = "#f0bc00", "TFi4" = "#5ac000")
my.pal.qc2 = c("FIB" = "#011736","Ci1" = "#400040", "Ci2" = "#87075c", "TFi1" = "#a11010", 
               "TFi2" = "#ad5100", "TFi3" = "#a38000", "TFi4" = "#367300")

#plot ratios
ggplot(dat.ratio.m, aes(Group, Ratio, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, outlier.size = 2, outlier.alpha = 1) +
  scale_fill_manual(values = my.pal.qc) +
  scale_color_manual(values = my.pal.qc2) +
  ylim(0, 60) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey70"),
        axis.ticks.x = element_blank())

#### FIGURE S4B ####

#calculate sum of counts mapped to genes
dat.sum = as.data.frame((colSums(dat.all[["tc.counts.qc"]]) + colSums(dat.all[["ercc.counts.qc"]]))/1e6) %>%
  mutate(Group = unlist(lapply(strsplit(row.names(.), "[.]"), function(X) X[[1]]))) %>%
  rownames_to_column(var="Cell")
names(dat.sum) = c("Cell", "Sum", "Group")

#melt data frame
dat.sum.m = melt(dat.sum, id=c("Sum", "Group")) %>%
  mutate(Group = factor(.$Group, levels = c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4")))

#plot sums
ggplot(dat.sum.m, aes(Group, Sum, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, outlier.size = 2, outlier.alpha = 1) +
  scale_fill_manual(values = my.pal.qc) +
  scale_color_manual(values = my.pal.qc2) +
  ylim(0, 7) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey70"),
        axis.ticks.x = element_blank())

#### FIGURE S4C ####

#calculate number of detected genes
dat.gene = data.frame(Cell = cells.qc, N.genes = colSums(dat.all[["tc.counts.qc"]] >= 1) / 1000, Group = dat.ratio$Group) %>%
  melt(id = c("N.genes", "Group")) %>%
  mutate(Group = factor(.$Group, levels = c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4")))

#plot number of genes
ggplot(dat.gene, aes(Group, N.genes, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, outlier.size = 2, outlier.alpha = 1) +
  scale_fill_manual(values = my.pal.qc) +
  scale_color_manual(values = my.pal.qc2) +
  ylim(0, 25) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey70"),
        axis.ticks.x = element_blank())

#### FIGURE S4D ####

samples = as.character(unique(dat.gene$Group))

#calculate CV for genes
dat.cv.tc = matrix(nrow = 0, ncol = 3)
for (i in 1:length(samples)) {
  
  dat.multi = dat.all[["tc.tpm.qc.log2"]][, grep(samples[i], colnames(dat.all[["tc.tpm.qc.log2"]]))] %>%
    .[which(rowMeans(.) > 0), ]
  dat.sd = apply(dat.multi, 1, sd)
  dat.mean = rowMeans(dat.multi)
  dat.cv = log10(dat.sd/dat.mean+1)
  dat.comb = as.data.frame(cbind(log10(dat.mean+1), dat.cv)) %>%
    mutate(Group = samples[i])
  dat.cv.tc = rbind(dat.cv.tc, dat.comb)}
dat.cv.tc[is.na(dat.cv.tc)] = 0

#calculate CV for ERCC spike-ins
dat.cv.ercc = matrix(nrow  =0, ncol = 3)
for (i in 1:length(samples)) {
  dat.ercc.multi = dat.all[["ercc.tpm.qc.log2"]][, grep(samples[i], colnames(dat.all[["ercc.tpm.qc.log2"]]))] %>%
    .[which(rowMeans(.) > 0), ]
  dat.ercc.sd = apply(dat.ercc.multi, 1, sd)
  dat.ercc.mean = rowMeans(dat.ercc.multi)
  dat.cv = log10(dat.ercc.sd/dat.ercc.mean+1)
  dat.comb.ercc = as.data.frame(cbind(log10(dat.ercc.mean+1), dat.cv)) %>%
    mutate(Group = samples[i])
  dat.cv.ercc = rbind(dat.cv.ercc, dat.comb.ercc)}
dat.cv.ercc[is.na(dat.cv.ercc)] = 0

#create color palette
my.pal.cv1 = c("FIB" = "#033782", "Ci1" = "#8b008b", "Ci2" = "#d40b91", "TFi1" = "#ec1717", 
               "TFi2" = "#fa7600", "TFi3" = "#f0bc00", "TFi4" = "#5ac000")
my.pal.cv2 = c("FIB" = "#011736", "Ci1" = "#400040", "Ci2" = "#87075c", "TFi1" = "#a11010", 
               "TFi2" = "#ad5100", "TFi3" = "#a38000", "TFi4" = "#367300")

#plot CV
ggplot(dat.cv.tc, aes(V1, dat.cv)) +
  geom_point(color = "grey80", alpha = 1) +
  geom_point(data = dat.cv.ercc, aes(V1, dat.cv, fill = Group, color = Group), pch = 21, size = 2, alpha = 1) +
  scale_fill_manual(values = my.pal.cv1) +
  scale_color_manual(values = my.pal.cv2) +
  geom_smooth(aes(V1, dat.cv, color = Group), se = F, size = 0.5, alpha = 0.5) +
  theme_jo()

#### FIGURE S4E ####

#read TPM data
dat.tc.ref = as.data.frame(as.matrix(tcs.iN@data))

#dendrogram 
dat.tpm = dat.tc.ref[which(row.names(dat.tc.ref) %in% tcs.iN@var.genes), which(colnames(dat.tc.ref) %in% tcs.iN.df$Cell)]

#create dendrogram
dat.dend = t(dat.tpm) %>%
  scale %>%
  dist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

#modify dendrogram
dend.label = dat.dend %>%
  labels %>%
  as.data.frame
dend.annot = data.frame(Cluster = as.character(tcs.iN.df$Cluster)) %>%
  mutate(rownames = tcs.iN.df$Cell,
         Group = tcs.iN.df$Sample,
         Cluster = chartr("12345", "ABCDE", Cluster))
dend.col = left_join(dend.annot, dend.label, by = c("rownames" = ".")) %>%
  mutate(Color2 = gsub2(c("A", "B", "C", "D", "E"), c("#033782", "#8b008b", "#d40b91", 
                                                      "#fa7600","#ec1717"), Cluster)) %>%
  column_to_rownames(var = "rownames")
dend.label = dend.label %>% column_to_rownames(var = ".")
dend.col = dend.col[row.names(dend.label), ] 
dend.ord = dend.col$Color

#plot dendrogram
dat.dend %>% 
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.8) %>%
  hang.dendrogram %>%
  set("leaves_col", dend.ord) %>%
  rotate(as.character(dend.annot[order(dend.annot$Cluster),]$rownames)) %>%
  set("labels", rep(NA, times=length(labels(.)))) %>%
  set("branches_lwd", 1.1) %>%
  set("branches_k_col", c("#033782", "#033782", "#8b008b", "#fa7600", "#d40b91", "#ec1717"), k=6) %>%
  hang.dendrogram(hang_height = 30) %>%
  plot

#with batch information
dend.label = dat.dend %>%
  labels %>%
  as.data.frame
dend.annot = data.frame(Cluster = as.character(tcs.iN.df$Cluster)) %>%
  mutate(rownames = tcs.iN.df$Cell,
         Cluster = chartr("12345", "ABCDE", Cluster),
         Group = tcs.iN.df$Sample)
dend.col = left_join(dend.annot, dend.label, by = c("rownames" = ".")) %>%
  mutate(Color = gsub2(c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4"), c("grey30", "grey30", "grey30", "#ec1717", "#fa7600", "#f0bc00", "#5ac000"), Group)) %>%
  column_to_rownames(var = "rownames")
dend.label = dend.label %>% column_to_rownames(var = ".")
dend.col = dend.col[row.names(dend.label), ] 
dend.ord = dend.col$Color

#plot
dat.dend %>% 
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.8) %>%
  hang.dendrogram %>%
  set("leaves_col", dend.ord) %>%
  rotate(as.character(dend.annot[order(dend.annot$Cluster),]$rownames)) %>%
  set("labels", rep(NA, times=length(labels(.)))) %>%
  set("branches_lwd", 1.1) %>%
  set("branches_col", "grey30") %>%
  hang.dendrogram(hang_height = 30) %>%
  plot