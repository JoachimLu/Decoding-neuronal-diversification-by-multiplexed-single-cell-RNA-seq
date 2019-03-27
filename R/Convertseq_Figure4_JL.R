#### FIGURE 4 ####

#load libraries
pacman::p_load(rio, tidyverse, data.table, Seurat, readxl, GOfuncR, pheatmap)

#custom functions
#--------------------------Read data------------------------------
readat <- function(id) {
  pid = paste0("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Figure.4/", id)
  dat = import(pid)
  dat = dat %>% remove_rownames %>% column_to_rownames(var = names(dat)[1])
}

readex <- function(id, sheet) {
  pid = paste0("C:/Users/joachim/Desktop/Projects/scRNAseq.iN/Figure.Data/Figure.4/", id)
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

darken <- function(color, factor=1.5){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) <- names(color)
  col
}

#### FIGURE 4B ####

#read.data
dat.ref.tc.count = readat("tc.mat.count.qc.csv")

#Seurat
tcs.iN = CreateSeuratObject(raw.data = dat.ref.tc.count,
                            min.cells = 3,
                            normalization.method = "LogNormalize") %>%
  FindVariableGenes(x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1) %>%
  ScaleData(genes.use = .@var.genes) %>%
  RunPCA() %>%
  ProjectPCA()%>%
  FindClusters(dims.use = 1:9, force.recalc = T) %>%
  RunTSNE(dims.use = 1:9, do.fast = T)
#TSNEPlot(tcs.iN, do.label = T)

current.cluster.ids = c(0:4)
new.cluster.ids = c("CL4", "CL5", "CL1", "CL2", "CL3")
tcs.iN@ident <- plyr::mapvalues(x = tcs.iN@ident, from = current.cluster.ids, to = new.cluster.ids)
tcs.iN@ident = factor(tcs.iN@ident, levels = c("CL1", "CL2", "CL3", "CL4", "CL5"))
tcs.iN@meta.data$res.0.8 = chartr("23401", "12345", tcs.iN@meta.data$res.0.8)

#create data frame for plotting
tcs.iN.df = as.data.frame(tcs.iN@dr$tsne@cell.embeddings) %>%
  rownames_to_column(var="Cell") %>%
  mutate(Sample = unlist(lapply(strsplit(Cell, "[.]"), function(X) X[[1]])),
         Cluster = factor(tcs.iN@meta.data$res.0.8, levels = c(1, 2, 3, 4, 5)))

#create color palettes
my.pal.c.fill = c("#033782", "#8b008b", "#d40b91", "#fa7600", "#ec1717")
my.pal.c.col = darken(my.pal.c.fill)

#plot
ggplot(tcs.iN.df, aes(tSNE_1, tSNE_2, fill = as.factor(Cluster), color = as.factor(Cluster), label = Cell)) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = my.pal.c.fill) +
  scale_color_manual(values = my.pal.c.col) +
  geom_text(size = 1) +
  theme_jo()

#### FIGURE 4C ####

#get cluster markers
cluster.markers = FindAllMarkers(tcs.iN, only.pos = T)
cluster.markers$pval.log = -log10(cluster.markers$p_val)

#create data frame for plotting heatmap
top10 = cluster.markers %>% group_by(cluster) %>% top_n(20, pval.log)
dat.heat = DoHeatmap(object = tcs.iN,
                     genes.use = top10$gene,
                     slim.col.label = TRUE,
                     remove.key = F)
dat.heat.df = dat.heat$data[, -4]
dat.heat.df = dcast(dat.heat.df, gene~cell) %>%
  column_to_rownames(var = "gene")
gene.ident = top10 %>%
  dplyr::select(c(gene, cluster, pval.log))
cell.ident = tcs.iN.df %>%
  dplyr::select(c(Cell, Cluster, Sample)) %>%
  mutate(Ident = as.numeric(as.character(Cluster))) %>%
  arrange(Ident)
m = data.frame(gene = row.names(dat.heat.df)) %>%
  left_join(gene.ident, by = "gene") %>%
  arrange(cluster, desc(pval.log))
dat.heat.df = dat.heat.df[m$gene, cell.ident$Cell]
pheatmap(as.matrix(dat.heat.df), cluster_cols = F, cluster_rows = F)

#create color annotation
dat.exo = cbind(data.frame(Sample = as.character(cell.ident$Sample)), 
                data.frame(Cluster = chartr("12345", "ABCDE", as.character(cell.ident$Cluster))))
row.names(dat.exo) = colnames(dat.heat.df)
ann_colors = list(
  Cluster = c("A" = "#033782", 
              "B" = "#8b008b", 
              "C" = "#d40b91", 
              "D" = "#fa7600", 
              "E" = "#ec1717"),
  Sample = c("FIB" = "#003829", 
             "Ci1" = "#385431", 
             "TFi1" = "#8b721e",
             "TFi2" = "#8b721e",
             "Ci2" = "#668414", 
             "TFi3" = "#ccc200",
             "TFi4" = "#ccc200"))

#plot heatmap
pheatmap(as.matrix(dat.heat.df), 
         border_color = "grey10",
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_col = dat.exo, 
         annotation_colors = ann_colors,
         annotation_legend = F,
         color = colorRampPalette(c("white", "#6d81a3", "#021734", "#fed38a", "#e99400"))(100),
         show_rownames = T,
         show_colnames = F,
         fontsize = 15,
         fontsize_row = 6,
         gaps_col = c(86, 171, 248, 359),
         gaps_row = c(16, 29, 38, 46))

#### FIGURE 4D ####

#read both datasets
cort = read.table(file = '/home/joachim/CORTEX/geneMatrix.tsv', sep = '\t', header = T)
dat.ref.tc.count = read.csv("/home/joachim/CORTEX/count.mat.tc.qc.csv", header = T, row.names = 1)
cort = cort %>%
  remove_rownames %>%
  column_to_rownames(var = "geneId")
dat.meta = read.csv("/home/joachim/CORTEX/cort.meta.csv", header = T, row.names = 1) %>%
  rownames_to_column(var = "ID")
id.cort = read.csv("/home/joachim/CORTEX/identifiers.csv", header = T, row.names = 1) %>%
  rownames_to_column(var = "ID") %>%
  .$ID %>%
  as.character
id.tc = colnames(dat.ref.tc.count)

#remove low expressed genes
dat.ref.tc.count = dat.ref.tc.count[which(rowSums(dat.ref.tc.count > 0) > 10), ]

#adjust genes and expression of both datasets
tc = cpm(dat.ref.tc.count)
#write.csv(tc, "tc.cpm.csv")
cort.id = cort[, which(colnames(cort) %in% id.cort)]
cort.id = cort.id[which(rowSums(cort.id > 0) > 10), ]
cort.id.raw = 2^cort.id
cort.id.raw[cort.id.raw == 1] = 0
cort.id.raw = cort.id.raw[which(row.names(cort.id.raw) %in% row.names(dat.ref.tc.count)), ]
tc = tc[which(row.names(tc) %in% row.names(cort.id.raw)), ]
cort.id.raw = cort.id.raw[row.names(tc), ]

#create Seurat objects 
dat.cort = CreateSeuratObject(raw.data = cort.id.raw, display.progress = T) %>%
  NormalizeData %>%
  ScaleData(display.progress = T) %>%
  FindVariableGenes

dat.tc = CreateSeuratObject(raw.data = tc, display.progress = T) %>%
  NormalizeData %>%
  ScaleData(display.progress = T) %>%
  FindVariableGenes

#take union of variable genes
hvg.cort <- rownames(x = head(x = dat.cort@hvg.info, n = 2000))
hvg.tc <- rownames(x = head(x = dat.tc@hvg.info, n = 2000))
hvg.union <- union(x = hvg.cort, y = hvg.tc)
dat.tc@meta.data[, "protocol"] <- "tc"
dat.cort@meta.data[, "protocol"] <- "cort"

#RunCCA
cort.tc <- RunCCA(object = dat.tc, object2 = dat.cort, genes.use = hvg.union)
cort.cor = as.matrix(cort.tc@data)

#exclude unwanted cells
exi = c("EN-PFC1", "EN-PFC2", "EN-PFC3", "EN-V1-1", "EN-V1-2", "EN-V1-3")
inhi = c("IN-CTX-CGE1", "IN-CTX-CGE2", "IN-CTX-MGE1", "IN-CTX-MGE2")
div = c("IPC-div1", "IPC-div2", "MGE-div", "RG-div1", "RG-div2")
prog = c("MGE-IPC1", "MGE-IPC2", "MGE-IPC3", "MGE-RG1", "MGE-RG2", "RG-early", "oRG", "vRG", "tRG")
endo = "Endothelial"
mural = "Mural"
microglia = "Microglia"
newborn = c("nIN1", "nIN2", "nIN3", "nIN4", "nIN5")
unknown = c("U1", "U2", "U3", "U4")
other = c("Astrocyte", "Choroid", "Glyc", "IPC-nEN1", "IPC-nEN2", "IPC-nEN3", "IN-STR", "nEN-early1", "nEN-early2", "nEN-late", "OPC")
dat.meta$Group = as.character(dat.meta$Class)
dat.meta$Group[which(dat.meta$Group %in% exi)] = "exi"
dat.meta$Group[which(dat.meta$Group %in% inhi)] = "inhi"
dat.meta$Group[which(dat.meta$Group %in% div)] = "div"
dat.meta$Group[which(dat.meta$Group %in% prog)] = "prog"
dat.meta$Group[which(dat.meta$Group %in% endo)] = "endo"
dat.meta$Group[which(dat.meta$Group %in% mural)] = "mural"
dat.meta$Group[which(dat.meta$Group %in% microglia)] = "microglia"
dat.meta$Group[which(dat.meta$Group %in% newborn)] = "newborn"
dat.meta$Group[which(dat.meta$Group %in% unknown)] = "unknown"
dat.meta$Group[which(dat.meta$Group %in% other)] = "other"
my.exi = dat.meta[which(dat.meta$Group == "exi"), ]$ID
my.inhi = dat.meta[which(dat.meta$Group == "inhi"), ]$ID
my.div = dat.meta[which(dat.meta$Group == "div"), ]$ID
my.prog = dat.meta[which(dat.meta$Group == "prog"), ]$ID
my.endo = dat.meta[which(dat.meta$Group == "endo"), ]$ID
my.mural = dat.meta[which(dat.meta$Group == "mural"), ]$ID
my.microglia = dat.meta[which(dat.meta$Group == "microglia"), ]$ID
my.newborn = dat.meta[which(dat.meta$Group == "newborn"), ]$ID
my.unknown = dat.meta[which(dat.meta$Group == "unknown"), ]$ID
my.other = dat.meta[which(dat.meta$Group == "other"), ]$ID
my.empty = dat.meta[which(dat.meta$Group == "empty"), ]$ID
tc.names = colnames(tc)
my.cells = c(my.exi, my.inhi, my.endo, my.mural, my.microglia, my.newborn, my.prog, my.div, tc.names)
cort.cor.2 = cort.cor[ , which(colnames(cort.cor) %in% my.cells)]

#correlate cells, keep cells with correlaiton >= 0.4
dat.cor = cor(cort.cor.2)
dat.cor.2 = dat.cor
diag(dat.cor.2) = 0
dat.cor.2[lower.tri(dat.cor.2)] = 0
dat.cor.2[dat.cor.2<0.3] = 0
dat.cor.2[dat.cor.2>=0.3] = 1

#compile correlation table for Cytoscape
my.cor = matrix(nrow = 0, ncol = 3)
for (i in 1:length(colnames(dat.cor.2))) {
  
  cor.names = row.names(dat.cor.2[which(dat.cor.2[,i] == 1), ])
  if (length(cor.names) > 0) {
    cor.names = data.frame("Target" = cor.names)
    cor.names$Cell = colnames(dat.cor.2)[i]
    cor.names$CORR = dat.cor[which(row.names(dat.cor) %in% cor.names$Target), i]
    colnames(cor.names) = c("Target", "Cell", "CORR")
    my.cor = rbind(my.cor, cor.names) }
  print(i)}
out1 = which(my.cor$Target %in% id.cort & my.cor$Cell %in% id.cort)
out2 = which(my.cor$Target %in% id.tc & my.cor$Cell %in% id.tc)
my.cor2 = my.cor
my.cor2 = my.cor2[-c(out2), ]
dat.meta = dat.meta %>%
  select(ID, Group)
my.cor2 = my.cor2 %>%
  left_join(., dat.meta, by = c("Target" = "ID"))
my.cor2 = my.cor2 %>%
  left_join(., dat.meta, by = c("Cell" = "ID"))
colnames(my.cor2) = c("Source", "Target", "CORR", "Group.Source", "Group.Target")

#### FIGURE 4E ####

#read data
dat.tc.ref = readat("tc.mat.tpm.qc.log2.csv")

#define exogenous genes
marker.genes.exo = c("ASCL1_EXO", "DLX1_EXO", "DLX2_EXO", "FEV_EXO", "FOXA2_EXO", "FOXP2_EXO", "ISL1_EXO", 
                     "LHX2_EXO", "NEUROD1_EXO", "NEUROG2_EXO", "NR2F1_EXO", "NR2F2_EXO", "NR4A2_EXO", 
                     "OLIG2_EXO", "OTX2_EXO", "PAX6_EXO", "PITX3_EXO", "POU3F2_EXO", "TLX3_EXO", "ZIC1_EXO")

#extract TPM of EXO marker genes early
final.dat.exo = matrix(ncol = 7, nrow = 0)

for (i in 1:length(marker.genes.exo)) {
  m.gene.pos = which(dat.tc.ref[marker.genes.exo[i], ] > 0)
  if(length(m.gene.pos) > 0) {
    dat.gene = tcs.iN.df[which(tcs.iN.df$Cell %in% colnames(dat.tc.ref[, m.gene.pos, drop = F])), ] %>%
      mutate(Gene = gsub("_EXO", "", marker.genes.exo[i]),
             TPM = as.numeric(dat.tc.ref[marker.genes.exo[i], which(colnames(dat.tc.ref) %in% .$Cell)]))
    final.dat.exo = rbind(final.dat.exo, dat.gene)}}

#scale
final.dat.exo = final.dat.exo %>%
  as.data.table %>%
  .[order(Gene, TPM)] %>%
  mutate(TPM = .[, scale(TPM), by = Gene]$V1)

#plot EXO data in facet grid
ggplot(tcs.iN.df, aes(tSNE_1, tSNE_2)) +
  geom_point(size = 2, color = "#081d3a") +
  geom_point(data = final.dat.exo, aes(tSNE_1, tSNE_2, color = TPM), size = 2) +
  scale_color_gradientn(colors = c("#081d3a", "#263b5a", "white","#e99400", "#ec1717")) +
  scale_x_continuous(breaks = c(-7, 5, 17)) +
  theme_jo() +
  theme(axis.text = element_blank(),
        strip.text = element_text(size = 22),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.text = element_text()) +
  facet_wrap(~Gene, nrow = 4)