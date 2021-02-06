#check if pacman package is loaded and install if not
if (!require(pacman)) {
  install.packages("pacman")
  if(!require(pacman)) stop("Package pacman not found")
}

#set data path
my_path <- "D:/Projects/Single_Cell_Convert_Seq"

#load libraries and custom functions
pacman::p_load(rio, tidyverse, pheatmap, data.table, Seurat, readxl, GOfuncR, Biobase, edgeR, monocle, dendextend)
source(file.path(my_path, "R", "utilities.R"))

#check package version of Seurat
if(package.version("Seurat") != "2.3.4") {
  stop("Analysis in this script is based on Seurat version 2.3.4 and wont run on newer versions")
}

if(package.version("monocle") != "2.2.0") {
  stop("Analysis in this script is based on monocle version 2.2.0 and wont run on newer versions")
}

####################################################################################################
############################################ FIGURE S4A ############################################
####################################################################################################

#read tc data
my_files <- c(file.path(my_path, "Data", "tc_mat_count_raw.csv"),
              file.path(my_path, "Data", "tc_mat_tpm_raw.csv"),
              file.path(my_path, "Data", "ercc_counts_raw.csv"),
              file.path(my_path, "Data", "ercc_tpm_raw.csv"),
              file.path(my_path, "Data", "tc_qc_low.csv"))
dat_all <- pbapply::pblapply(my_files, readat)
names(dat_all) <- basename(gsub("\\.csv$", "", my_files))

#remove dead cells and cell doublets from count data
dat_count <- dat_all[["tc_mat_count_raw"]][ ,-(which(colnames(dat_all[["tc_mat_count_raw"]]) %in% dat_all[["tc_qc_low"]]$LowQC))]

#remove cells with low number of detected genes
gene_count <- apply(dat_count, 2, function(y) length(which(y > 0)))
dat_all[["tc_counts_qc"]] <- as.data.frame(dat_count[, -(which(gene_count < (mean(gene_count) - 3*sd(gene_count))))])
cells_qc <- names(dat_all[["tc_counts_qc"]])
dat_all[["tc_tpm_qc"]] <- dat_all[["tc_mat_tpm_raw"]][, which(names(dat_all[["tc_mat_tpm_raw"]]) %in% cells_qc)]
dat_all[["ercc_counts_qc"]] <- dat_all[["ercc_counts_raw"]][, which(names(dat_all[["ercc_counts_raw"]]) %in% cells_qc)]
dat_all[["ercc_tpm_qc"]] <- dat_all[["ercc_tpm_raw"]][, which(names(dat_all[["ercc_tpm_raw"]]) %in% cells_qc)]
dat_all[["tc_counts_qc_log2"]] <- log2(dat_all[["tc_counts_qc"]] + 1)
dat_all[["tc_tpm_qc_log2"]] <- log2(dat_all[["tc_tpm_qc"]] + 1)
dat_all[["ercc_tpm_qc_log2"]] = log2(dat_all[["ercc_tpm_qc"]] + 1)

#calculate percentage of ERCC spikes
dat_ratio <- as.data.frame((100*(colSums(dat_all[["ercc_counts_qc"]])) / 
                             (colSums(dat_all[["ercc_counts_qc"]]) + colSums(dat_all[["tc_counts_qc"]])))) %>%
  mutate(Group <- unlist(lapply(strsplit(row.names(.), "[.]"), function(X) X[[1]]))) %>%
  rownames_to_column(var = "Cell")
names(dat_ratio) <- c("Cell", "Ratio", "Group")

#melt data frame
dat_ratio_m <- melt(dat_ratio, id = c("Ratio", "Group")) %>%
  mutate(Group = factor(.$Group, levels = c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4")))

#create color palette
my_pal_qc <- c("FIB" = "#033782", "Ci1" = "#8b008b", "Ci2" = "#d40b91", "TFi1" = "#ec1717", 
              "TFi2" = "#fa7600", "TFi3" = "#f0bc00", "TFi4" = "#5ac000")
my_pal_qc2 <- darken(my_pal_qc)

#plot ratios
ggplot(dat_ratio_m, aes(Group, Ratio, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, outlier.size = 2, outlier.alpha = 1) +
  scale_fill_manual(values = my_pal_qc) +
  scale_color_manual(values = my_pal_qc2) +
  ylim(0, 60) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey70"),
        axis.ticks.x = element_blank())

####################################################################################################
############################################ FIGURE S4B ############################################
####################################################################################################

#calculate sum of counts mapped to genes
dat_sum <- as.data.frame((colSums(dat_all[["tc_counts_qc"]]) + colSums(dat_all[["ercc_counts_qc"]]))/1e6) %>%
  mutate(Group = unlist(lapply(strsplit(row.names(.), "[.]"), function(X) X[[1]]))) %>%
  rownames_to_column(var = "Cell")
names(dat_sum) <- c("Cell", "Sum", "Group")

#melt data frame
dat_sum_m <- melt(dat_sum, id = c("Sum", "Group")) %>%
  mutate(Group = factor(.$Group, levels = c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4")))

#plot sums
ggplot(dat_sum_m, aes(Group, Sum, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, outlier.size = 2, outlier.alpha = 1) +
  scale_fill_manual(values = my_pal_qc) +
  scale_color_manual(values = my_pal_qc2) +
  ylim(0, 7) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey70"),
        axis.ticks.x = element_blank())

####################################################################################################
############################################ FIGURE S4C ############################################
####################################################################################################

#calculate number of detected genes
dat_gene <- data.frame(Cell = cells_qc, 
                       N.genes = colSums(dat_all[["tc_counts_qc"]] >= 1) / 1000, 
                       Group = dat_ratio$Group) %>%
  melt(id = c("N.genes", "Group")) %>%
  mutate(Group = factor(.$Group, levels = c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4")))

#plot number of genes
ggplot(dat_gene, aes(Group, N.genes, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, outlier.size = 2, outlier.alpha = 1) +
  scale_fill_manual(values = my_pal_qc) +
  scale_color_manual(values = my_pal_qc2) +
  ylim(0, 25) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey70"),
        axis.ticks.x = element_blank())

####################################################################################################
############################################ FIGURE S4D ############################################
####################################################################################################

samples <- as.character(unique(dat_gene$Group))

#calculate CV for genes
dat_cv_tc <- matrix(nrow = 0, ncol = 3)
for (i in 1:length(samples)) {
  
  dat_multi <- dat_all[["tc_tpm_qc_log2"]][, grep(samples[i], colnames(dat_all[["tc_tpm_qc_log2"]]))] %>%
    .[which(rowMeans(.) > 0), ]
  dat_sd <- apply(dat_multi, 1, sd)
  dat_mean <- rowMeans(dat_multi)
  dat_cv <- log10(dat_sd/dat_mean + 1)
  dat_comb <- as.data.frame(cbind(log10(dat_mean+1), dat_cv)) %>%
    mutate(Group = samples[i])
  dat_cv_tc <- rbind(dat_cv_tc, dat_comb)}
dat_cv_tc[is.na(dat_cv_tc)] = 0

#calculate CV for ERCC spike-ins
dat_cv_ercc <- matrix(nrow = 0, ncol = 3)
for (i in 1:length(samples)) {
  dat_ercc_multi <- dat_all[["ercc_tpm_qc_log2"]][, grep(samples[i], colnames(dat_all[["ercc_tpm_qc_log2"]]))] %>%
    .[which(rowMeans(.) > 0), ]
  dat_ercc_sd <- apply(dat_ercc_multi, 1, sd)
  dat_ercc_mean <- rowMeans(dat_ercc_multi)
  dat_cv <- log10(dat_ercc_sd/dat_ercc_mean + 1)
  dat_comb_ercc <- as.data.frame(cbind(log10(dat_ercc_mean + 1), dat_cv)) %>%
    mutate(Group = samples[i])
  dat_cv_ercc <- rbind(dat_cv_ercc, dat_comb_ercc)}
dat_cv_ercc[is.na(dat_cv_ercc)] = 0

#create color palette
my_pal_cv1 <- c("FIB" = "#033782", "Ci1" = "#8b008b", "Ci2" = "#d40b91", "TFi1" = "#ec1717", 
                "TFi2" = "#fa7600", "TFi3" = "#f0bc00", "TFi4" = "#5ac000")
my_pal_cv2 <- c("FIB" = "#011736", "Ci1" = "#400040", "Ci2" = "#87075c", "TFi1" = "#a11010", 
                "TFi2" = "#ad5100", "TFi3" = "#a38000", "TFi4" = "#367300")

#plot CV
ggplot(dat_cv_tc, aes(V1, dat_cv)) +
  geom_point(color = "grey80", alpha = 1) +
  geom_point(data = dat_cv_ercc, aes(V1, dat_cv, fill = Group, color = Group), pch = 21, size = 2, alpha = 1) +
  scale_fill_manual(values = my_pal_cv1) +
  scale_color_manual(values = my_pal_cv2) +
  geom_smooth(aes(V1, dat_cv, color = Group), se = F, size = 0.5, alpha = 0.5) +
  theme_jo()

####################################################################################################
############################################ FIGURE S4E ############################################
####################################################################################################

#read count 
dat_ref_tc_count <- readat(file.path(my_path, "Data", "tc_mat_count_qc.csv"))

#Seurat
tcs_iN <- CreateSeuratObject(raw.data = dat_ref_tc_count,
                             min.cells = 3,
                             normalization.method = "LogNormalize") %>%
  FindVariableGenes(x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1) %>%
  ScaleData(genes.use = .@var.genes) %>%
  RunPCA() %>%
  ProjectPCA()%>%
  FindClusters(dims.use = 1:9, force.recalc = T) %>%
  RunUMAP(dims.use = 1:10)

#create data frame for ggplot
tcs_iN_df <- as.data.frame(tcs_iN@dr$umap@cell.embeddings) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Cluster = factor(tcs_iN@meta.data$res.0.8),
         Sample = strsplit2(Cell, "[.]")[,1])

my_pal_fill <- c("#3b322a", "#8f856e", "#ffedc5", "#ffb73d", "#ef5134")
my_pal_col <- darken(my_pal_fill)

#plot UMAP with cluster information
ggplot(tcs_iN_df, aes(UMAP1, UMAP2, fill = as.factor(Cluster), color = as.factor(Cluster))) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = my_pal_fill) +
  scale_color_manual(values = my_pal_col) +
  theme_jo()

#extract reference data
dat_tc_ref <- as.data.frame(as.matrix(tcs_iN@data))

#dendrogram 
dat_tpm <- dat_tc_ref[which(row.names(dat_tc_ref) %in% tcs_iN@var.genes), which(colnames(dat_tc_ref) %in% tcs_iN_df$Cell)]

#create dendrogram
dat_dend <- t(dat_tpm) %>%
  scale %>%
  dist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

#modify dendrogram
dend_label <- dat_dend %>%
  labels %>%
  as.data.frame
dend_annot <- data.frame(Cluster = as.character(tcs_iN_df$Cluster)) %>%
  mutate(rownames = tcs_iN_df$Cell,
         Group = tcs_iN_df$Sample,
         Cluster = chartr("01234", "DEABC", Cluster))
dend_col <- left_join(dend_annot, dend_label, by = c("rownames" = ".")) %>%
  mutate(Color2 = gsub2(c("A", "B", "C", "D", "E"), c("#3b322a", "#8f856e", "#ffedc5", 
                                                      "#ffb73d","#ef5134"), Cluster)) %>%
  column_to_rownames(var = "rownames")
dend_label <- dend_label %>% column_to_rownames(var = ".")
dend_col <- dend_col[row.names(dend_label), ] 
dend_ord <- dend_col$Color

#plot dendrogram
dat_dend %>% 
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.8) %>%
  hang.dendrogram %>%
  set("leaves_col", dend_ord) %>%
  rotate(as.character(dend_annot[order(dend_annot$Cluster),]$rownames)) %>%
  set("labels", rep(NA, times = length(labels(.)))) %>%
  set("branches_lwd", 1.1) %>%
  set("branches_k_col", c("#3b322a", "#3b322a", "#8f856e",  "#ffb73d", "#ffedc5", "#ef5134"), k = 6) %>%
  hang.dendrogram(hang_height = 30) %>%
  plot

#with batch information
dend_label <- dat_dend %>%
  labels %>%
  as.data.frame
dend_annot <- data.frame(Cluster = as.character(tcs_iN_df$Cluster)) %>%
  mutate(rownames = tcs_iN_df$Cell,
         Group = tcs_iN_df$Sample,
         Cluster = chartr("01234", "DEABC", Cluster))
dend_col <- left_join(dend_annot, dend_label, by = c("rownames" = ".")) %>%
  mutate(Color = gsub2(c("FIB", "Ci1", "Ci2", "TFi1", "TFi2", "TFi3", "TFi4"), c("grey30", "grey30", "grey30", "#ec1717", "#fa7600", "#f0bc00", "#5ac000"), Group)) %>%
  column_to_rownames(var = "rownames")
dend_label <- dend_label %>% column_to_rownames(var = ".")
dend_col <- dend_col[row.names(dend_label), ] 
dend_ord <- dend_col$Color

#plot
dat_dend %>% 
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.8) %>%
  hang.dendrogram %>%
  set("leaves_col", dend_ord) %>%
  rotate(as.character(dend_annot[order(dend_annot$Cluster),]$rownames)) %>%
  set("labels", rep(NA, times = length(labels(.)))) %>%
  set("branches_lwd", 1.1) %>%
  set("branches_col", "grey30") %>%
  hang.dendrogram(hang_height = 30) %>%
  plot

####################################################################################################
############################################ FIGURE S4F ############################################
####################################################################################################

#plot UMAP with sample information
my_pal_fill <- c("FIB" = "grey30",
                 "Ci1" = "grey30", 
                 "Ci2" = "grey30", 
                 "TFi1" = "#ec1616",
                 "TFi2" ="#fa7602",
                 "TFi3" = "#f1bc00",
                 "TFi4" = "#5bc000")
my_pal_col <- darken(my_pal_fill)

ggplot(tcs_iN_df, aes(UMAP1, UMAP2, fill = as.factor(Sample), color = as.factor(Sample))) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = my_pal_fill) +
  scale_color_manual(values = my_pal_col) +
  theme_jo()

####################################################################################################
############################################ FIGURE S4G ############################################
####################################################################################################

iNM_sample_sheet <- readat(file.path(my_path, "Data", "tc_sample_sheet.csv"))
iNM_sample_sheet$Cluster <- tcs_iN_df$Cluster
iNM_gene_annotation <- readat(file.path(my_path, "Data", "tc_gene_annotation.csv"))
iNM_expr_matrix <- readat(file.path(my_path, "Data", "tc_mat_tpm_qc.csv"))
row.names(iNM_gene_annotation) <- row.names(iNM_expr_matrix)

#create cell data set
pd <- new("AnnotatedDataFrame", data = iNM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = iNM_gene_annotation)
iNM <- newCellDataSet(as.matrix(iNM_expr_matrix), phenoData = pd, featureData = fd)

#filter Genes
iNM <- detectGenes(iNM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(iNM), num_cells_expressed >= 3))

#diff_test <- differentialGeneTest(iNM[expressed_genes, ], fullModelFormulaStr = "~Media", cores = 4)
#sig_genes <- subset(diff_test, qval < 0.01) %>% filter(gene_short_name %in% expressed_genes)
#ordering_genes <- sig_genes$gene_short_name

#read differentially expressed gene data
ordering_genes <- readat(file.path(my_path, "Data", "tc_pt_unsup_diff_genes.csv")) %>% 
  .$gene_short_name

#reduce dimension, order cells
iNM <- setOrderingFilter(iNM, ordering_genes) %>%
  reduceDimension(reduction_method = "DDRTree") %>%
  orderCells(reverse = F)

#obtain data frame
pt_iN <- plot_cell_trajectory(iNM, color_by = "State", show_tree = T, show_branch_points = T)
pt_dat <- pt_iN$data %>%
  column_to_rownames(var = "sample_name") %>%
  .[row.names(as.data.frame(t(iNM@reducedDimS))), ]

#extract information for transferring data to ggplot2
iNM_coord <- as.data.frame(t(iNM@reducedDimS)) %>%
  mutate(Batch = unlist(lapply(strsplit(row.names(.), "[.]"), function(X) X[[1]])),
         Cell = row.names(as.data.frame(t(iNM@reducedDimS))),
         State = pt_dat$State,
         Branch = chartr("123456789", "AAACBBBBB", State),
         Sample = gsub("TFi2", "TFi1", Batch),
         Sample = gsub("TFi3", "TFi2", Sample),
         Sample = gsub("TFi4", "TFi2", Sample)) %>%
  left_join(., pt_iN$data[ ,c("sample_name", "Pseudotime")], by = c("Cell" = "sample_name"))

#obtain minimum spanning tree
reduced_dim_coords <- reducedDimK(iNM)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(1, 2), ])) %>%
  select(prin_graph_dim_1 = X1, prin_graph_dim_2 = X2) %>%
  rownames_to_column(var = "sample_name") %>%
  mutate(sample_state = sample_name)
iN_mst <- minSpanningTree(iNM)
edge_list <- as.data.frame(igraph::get.edgelist(iN_mst)) %>%
  select(source = V1, target = V2)
edge_df_iN <- full_join(ica_space_df, edge_list, by = c("sample_name" = "source")) %>%
  plyr::rename(c("prin_graph_dim_1" = "source_prin_graph_dim_1", "prin_graph_dim_2" = "source_prin_graph_dim_2")) %>%
  full_join(ica_space_df[ ,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by = c("target" = "sample_name")) %>%
  plyr::rename(c("prin_graph_dim_1" = "target_prin_graph_dim_1", "prin_graph_dim_2" = "target_prin_graph_dim_2"))

#create color palettes
my_pal_tc_monoc1 <- c("FIB" = "#033782", "Ci1" = "#8b008b", "Ci2" = "#c71585", 
                      "TFi1" = "#ee2c2c", "TFi2" = "#ee2c2c", "TFi3" = "#ee7600", "TFi4" = "#ee7600")
my_pal_tc_monoc2 <- c("FIB" = "#011736","Ci1" = "#400040", "Ci2" = "#7a0d52", 
                      "TFi1" = "#a11d1d", "TFi2" = "#a11d1d", "TFi3" = "#a15000", "TFi4" = "#a15000")

#plot with sample information
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"), size = 1.2, linetype = "solid", na.rm = TRUE, data = edge_df_iN) +
  geom_point(data = iNM_coord[sample(nrow(iNM_coord)), ], aes(V1, V2, fill = Batch, color = Batch), size = 7, shape = 21) +
  scale_fill_manual(values = my_pal_tc_monoc1) +
  scale_color_manual(values = my_pal_tc_monoc2) +
  theme_jo() +
  theme(panel.border = element_rect(fill = NA, size = 1.5, color = "grey10"),
        axis.line = element_blank())

#plot pseudotime information
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = 1.2, linetype = "solid", na.rm = TRUE, 
               data = edge_df_iN, color="grey20") +
  geom_point(data = iNM_coord, aes(V1, V2, fill = Pseudotime), size = 9, shape = 21) +
  scale_fill_gradient2(low = "#60b4b4", mid = "#004d81", high = "#f69c0d", midpoint = max(iNM_coord$Pseudotime/2)) +
  theme_blank() +
  theme(panel.border = element_rect(fill = NA, size = 3, color = "grey10"),
        axis.ticks = element_blank(),
        axis.line = element_blank())

#plot facet
ggplot(iNM_coord, aes(V1, V2)) +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"), size = 1.2, linetype = "solid", na.rm = TRUE, data = edge_df_iN) +
  stat_density_2d(aes(color = Sample)) +
  geom_point(aes(fill = Sample, color = Sample), size = 6, shape = 21) +
  scale_fill_manual(values = my_pal_tc_monoc1) +
  scale_color_manual(values  = my_pal_tc_monoc2) +
  theme_blank() +
  theme(panel.border = element_rect(fill = NA, size = 1.5, color = "grey10"),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~Sample)

####################################################################################################
############################################ FIGURE S4H ############################################
####################################################################################################

#use protein coding genes that are expressed in at least 20 cells
expressed_genes <- row.names(subset(fData(iNM), num_cells_expressed >= 20))
genes <- readat(file.path(my_path, "Data", "gene_class.csv")) %>%
  filter(Class == "protein_coding") %>%
  .$Gene %>%
  .[which(. %in% expressed_genes)]

#Pseudotime-dependent genes
#differential gene expression along pseudotime
diff_test_res_pt <- differentialGeneTest(iNM[genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 4)
sig_gene_names <- row.names(subset(diff_test_res_pt, qval < 1e-4))

#scale matrix
cds_subset <- iNM[sig_gene_names,]
newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                       max(pData(cds_subset)$Pseudotime), length.out = 100))
m <- genSmoothCurves(cds_subset, 
                     cores = 1, 
                     trend_formula = "~sm.ns(Pseudotime, df = 3)",  
                     relative_expr = T, 
                     new_data = newdata)
m <- m[!apply(m, 1, sum) == 0,]
m <- log10(m + 1)
m <- m[!apply(m, 1, sd) == 0, ]
m <- Matrix::t(scale(Matrix::t(m), center = TRUE))
m <- m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] <- 0
m[m > 3] <- 3
m[m < -3] <- -3
pt_matrix <- m

#create pheatmap object
ph <- pheatmap(pt_matrix, cluster_cols = F, cutree_rows = 3, clustering_method = "ward.D2", fontsize_row = 1)

#create color annotation
annot_row <- as.data.frame(cutree(ph$tree_row, k = 3)) %>%
  rownames_to_column(var = "Gene")
colnames(annot_row) = c("Gene", "Cluster")
annot_row$Cluster = chartr("123", "ACB", annot_row$Cluster)
dat_dend = as.dendrogram(ph$tree_row) %>%
  rotate(annot_row[order(annot_row$Cluster), ]$Gene)
annot_row = annot_row %>%
  column_to_rownames(var = "Gene")
my_clust = as.hclust(dat_dend)

ann_colors = list(Cluster = c("A" = "#5aaab1", "B" = "#1a5b87", "C" = "#fda24a"))

#plot heatmap
pheatmap(pt_matrix, cluster_rows = my_clust, cluster_cols = F, annotation_row = annot_row, annotation_colors = ann_colors, 
         show_rownames = T, show_colnames = F, annotation_legend = F, border_color = NA, fontsize_row = 0.5, cutree_rows = 3,
         fontsize = 20, clustering_method = "ward.D2", legend = F)

####################################################################################################
########################################### FIGURE S4J-L ###########################################
####################################################################################################

#read log-transformed tpm data
dat_tc_ref = readat(file.path(my_path, "Data", "tc_mat_tpm_qc_log2.csv"))

#plot genes in pseudotime
#create data frame
dat_pt <- iNM_coord %>%
  .[order(.$Pseudotime), ]

#set marker genes
m_genes_fib <- c("CCNB1", "MKI67", "TK1", "TOP2A")
m_genes_neuro <- c("NRCAM", "SFRP1", "SNAP25", "SYT1")
m_genes_blood <- c("VEGFA", "FAT4", "PGF", "BMP4")

#create data frames for plotting
#-------------------------------Fibroblast genes-----------------------------------
#create in data frame for plotting
dat_pt_final_fib <- matrix(nrow = 0, ncol = 4)

for (i in seq_len(length(m_genes_fib))) {
  
  dat_marker_fib <- as.data.frame(t(dat_tc_ref[m_genes_fib[i], ])) %>%
    mutate(Sample = iNM_coord$Sample,
           Pseudotime = iNM_coord$Pseudotime,
           Gene = m_genes_fib[i])
  names(dat_marker_fib) <- c("Tpm", "Sample", "Pseudotime", "Gene")
  dat_pt_final_fib <- rbind(dat_pt_final_fib, dat_marker_fib)
}

#-------------------------------Neuronal genes-----------------------------------
#create in data frame for plotting
dat_pt_final_neuro <- matrix(nrow = 0, ncol = 4)

for (i in seq_len(length(m_genes_neuro))) {
  
  dat_marker_neuro <- as.data.frame(t(dat_tc_ref[m_genes_neuro[i], ])) %>%
    mutate(Sample = iNM_coord$Sample,
           Pseudotime = iNM_coord$Pseudotime,
           Gene = m_genes_neuro[i])
  names(dat_marker_neuro) <- c("Tpm", "Sample", "Pseudotime", "Gene")
  dat_pt_final_neuro <- rbind(dat_pt_final_neuro, dat_marker_neuro)
}

#-------------------------------Vascular genes-----------------------------------
#create in data frame for plotting
dat_pt_final_blood = matrix(nrow = 0, ncol = 4)

for (i in seq_len(length(m_genes_blood))) {
  
  dat_marker_blood <- as.data.frame(t(dat_tc_ref[m_genes_blood[i], ])) %>%
    mutate(Sample = iNM_coord$Sample,
           Pseudotime = iNM_coord$Pseudotime,
           Gene = m_genes_blood[i])
  names(dat_marker_blood) <- c("Tpm", "Sample", "Pseudotime", "Gene")
  dat_pt_final_blood <- rbind(dat_pt_final_blood, dat_marker_blood)
}

#modify marker levels
dat_pt_final_fib$Gene <- factor(dat_pt_final_fib$Gene, levels = c("CCNB1", "MKI67", "TK1", "TOP2A"))

#create color palette
mypal_pt1 <- c("FIB" = "#033782", "Ci1" = "#8b008b", "Ci2" = "#c71585", "TFi1" = "#ee2c2c", "TFi2" = "#ee7600")
mypal_pt2 <- c("FIB" = "#011736","Ci1" = "#400040", "Ci2" = "#7a0d52", "TFi1" = "#a11d1d", "TFi2" = "#a15000")

#----------------------------FigureS4J--------------------------------
#plot
ggplot(dat_pt_final_fib, aes(Pseudotime, Tpm)) +
  geom_point(aes(color = as.factor(Sample), fill = as.factor(Sample)), size = 5, alpha = 1, shape = 21) +
  scale_color_manual(values = mypal_pt2) +
  scale_fill_manual(values = mypal_pt1) +
  geom_smooth(se = T, size = 2, color = "grey30", span = 0.5) +
  theme_jo()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  facet_wrap(~Gene, ncol = 1, scales = "free")

#----------------------------FigureS4K--------------------------------
#plot
ggplot(dat_pt_final_neuro, aes(Pseudotime, Tpm)) +
  geom_point(aes(color = as.factor(Sample), fill = as.factor(Sample)), size = 5, alpha = 1, shape = 21) +
  scale_color_manual(values = mypal_pt2) +
  scale_fill_manual(values = mypal_pt1) +
  geom_smooth(se = T, size = 2, color = "grey30", span = 0.5) +
  theme_jo()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  facet_wrap(~Gene, ncol = 1, scales = "free")

#----------------------------FigureS4L--------------------------------
#plot
ggplot(dat_pt_final_blood, aes(Pseudotime, Tpm)) +
  geom_point(aes(color = as.factor(Sample), fill = as.factor(Sample)), size = 5, alpha = 1, shape = 21) +
  scale_color_manual(values = mypal_pt2) +
  scale_fill_manual(values = mypal_pt1) +
  geom_smooth(se = T, size = 2, color = "grey30", span = 0.5) +
  theme_jo()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  facet_wrap(~Gene, ncol = 1, scales = "free")

####################################################################################################
############################################ FIGURE S4M ############################################
####################################################################################################

#run monocle on developmental genes
iNM_sample_sheet <- readat(file.path(my_path, "Data", "tc_sample_sheet.csv"))
iNM_sample_sheet$Cluster <- tcs_iN_df$Cluster
iNM_gene_annotation <- readat(file.path(my_path, "Data", "tc_gene_annotation.csv"))
iNM_expr_matrix <- readat(file.path(my_path, "Data", "tc_mat_tpm_qc.csv"))
row.names(iNM_gene_annotation) <- row.names(iNM_expr_matrix)

#create cell data set
pd <- new("AnnotatedDataFrame", data = iNM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = iNM_gene_annotation)
iNM <- newCellDataSet(as.matrix(iNM_expr_matrix), phenoData = pd, featureData = fd)

#filter Genes
iNM <- detectGenes(iNM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(iNM), num_cells_expressed >= 3))

#obtain nervous and blood ordering genes
ordering_genes <- readat(file.path(my_path, "Data", "go_development_all.csv")) %>%
  .$Gene
ordering_genes <- unique(ordering_genes)
ordering_genes <- expressed_genes[which(expressed_genes %in% ordering_genes)]

#reduce dimension, order cells
iNM <- setOrderingFilter(iNM, ordering_genes) %>%
  reduceDimension(reduction_method = "DDRTree") %>%
  orderCells(reverse = F)

#obtain data frame
pt_iN <- plot_cell_trajectory(iNM, color_by = "State", show_tree = T, show_branch_points = T)
pt_dat <- pt_iN$data %>%
  column_to_rownames(var = "sample_name") %>%
  .[row.names(as.data.frame(t(iNM@reducedDimS))), ]

#extract information for transferring data to ggplot2
iNM_coord <- as.data.frame(t(iNM@reducedDimS)) %>%
  mutate(Sample = unlist(lapply(strsplit(row.names(.), "[.]"), function(X) X[[1]])),
         Cell = row.names(as.data.frame(t(iNM@reducedDimS))),
         State = pt_dat$State,
         Branch = chartr("123456789", "AAACBBBBB", State)) %>%
  left_join(., pt_iN$data[ ,c("sample_name", "Pseudotime")], by = c("Cell" = "sample_name"))

#obtain minimum spanning tree
reduced_dim_coords <- reducedDimK(iNM)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(1, 2), ])) %>%
  select(prin_graph_dim_1 = X1, prin_graph_dim_2 = X2) %>%
  rownames_to_column(var = "sample_name") %>%
  mutate(sample_state = sample_name)
iN_mst <- minSpanningTree(iNM)
edge_list <- as.data.frame(igraph::get.edgelist(iN_mst)) %>%
  select(source = V1, target = V2)
edge_df_iN <- full_join(ica_space_df, edge_list, by = c("sample_name" = "source")) %>%
  plyr::rename(c("prin_graph_dim_1" = "source_prin_graph_dim_1", "prin_graph_dim_2" = "source_prin_graph_dim_2")) %>%
  full_join(ica_space_df[ ,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by = c("target" = "sample_name")) %>%
  plyr::rename(c("prin_graph_dim_1" = "target_prin_graph_dim_1", "prin_graph_dim_2" = "target_prin_graph_dim_2"))

#define exogenous genes
marker_genes_exo <- c("ASCL1_EXO", "DLX1_EXO", "DLX2_EXO", "FEV_EXO", "FOXA2_EXO", "FOXP2_EXO", "ISL1_EXO",
                      "LHX2_EXO", "NEUROD1_EXO", "NEUROG2_EXO", "NR2F1_EXO", "NR2F2_EXO", "NR4A2_EXO", 
                      "OLIG2_EXO", "OTX2_EXO", "PAX6_EXO", "PITX3_EXO", "POU3F2_EXO", "TLX3_EXO", "ZIC1_EXO")

#extract TPM of EXO marker genes early
final_dat_exo <- matrix(ncol = 7, nrow = 0)

for (i in 1:length(marker_genes_exo)) {
  m_gene_pos <- which(dat_tc_ref[marker_genes_exo[i], ] > 0)
  if(length(m_gene_pos) > 0) {
    dat_gene <- iNM_coord[which(iNM_coord$Cell %in% colnames(dat_tc_ref[, m_gene_pos, drop = F])), ] %>%
      mutate(Gene = gsub("_EXO", "", marker_genes_exo[i]),
             TPM = as.numeric(dat_tc_ref[marker_genes_exo[i], which(colnames(dat_tc_ref) %in% .$Cell)]))
    final_dat_exo = rbind(final_dat_exo, dat_gene)}}

#scale
final_dat_exo <- final_dat_exo %>%
  as.data.table %>%
  .[order(Gene, TPM)] %>%
  mutate(TPM = .[, scale(TPM), by = Gene]$V1)

#plot EXO data in facet grid early
ggplot(iNM_coord, aes(V1, V2)) +
  geom_point(size = 3, color = "#5a5a5a") +
  geom_point(data = final_dat_exo, aes(V1, V2, color = TPM), size = 3) +
  scale_color_gradientn(colors = c("#5a5a5a", "white", "#f2a247", "#ee2c2c")) +
  scale_x_continuous(breaks = c(-7, 5, 17)) +
  theme_jo() +
  theme(axis.text = element_blank(),
        strip.text = element_text(size = 22),
        panel.border = element_rect(fill = NA, size = 1.5, color = "grey10"),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~Gene, nrow = 3)