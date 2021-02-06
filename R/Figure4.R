#check if pacman package is loaded and install if not
if (!require(pacman)) {
  install.packages("pacman")
  if(!require(pacman)) stop("Package pacman not found")
}

#set data path
my_path <- "D:/Projects/Single_Cell_Convert_Seq"

#load libraries and custom functions
pacman::p_load(rio, tidyverse, pheatmap, data.table, Seurat, readxl, GOfuncR, Biobase, edgeR, monocle)
source(file.path(my_path, "R", "utilities.R"))

#check package version of Seurat
if(package.version("Seurat") != "2.3.4") {
  stop("Analysis in this script is based on Seurat version 2.3.4 and wont run on newer versions")
}

if(package.version("monocle") != "2.2.0") {
  stop("Analysis in this script is based on monocle version 2.2.0 and wont run on newer versions")
}

####################################################################################################
############################################ FIGURE 4B #############################################
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
  FindClusters(dims.use = 1:10, force.recalc = T) %>%
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

#plot UMAP with sample information
my_pal_fill <- c("FIB" = "#233d7d",
                "Ci1" = "#862887", 
                "Ci2" = "#c61783", 
                "TFi1" = "#ed2d2d",
                "TFi2" ="#ed2d2d",
                "TFi3" = "#ec7723",
                "TFi4" = "#ec7723")
my_pal_col <- darken(my_pal_fill)

ggplot(tcs_iN_df, aes(UMAP1, UMAP2, fill = as.factor(Sample), color = as.factor(Sample))) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = my_pal_fill) +
  scale_color_manual(values = my_pal_col) +
  theme_jo()

####################################################################################################
############################################ FIGURE 4C #############################################
####################################################################################################

#define exogenous genes
marker_genes_exo <- c("ASCL1_EXO","DLX1_EXO","DLX2_EXO","FEV_EXO","FOXA2_EXO", "FOXP2_EXO",
                      "ISL1_EXO", "LHX2_EXO", "NEUROD1_EXO", "NEUROG2_EXO", "NR2F1_EXO","NR2F2_EXO","NR4A2_EXO",
                      "OLIG2_EXO", "OTX2_EXO", "PAX6_EXO","PITX3_EXO","POU3F2_EXO", "TLX3_EXO", "ZIC1_EXO")

#extract TPM of EXO marker genes
final_dat_exo <- matrix(ncol = 7, nrow = 0)
for (i in 1:length(marker_genes_exo)) {
  m_gene_pos <- which(dat_ref_tc_count[marker_genes_exo[i], ] > 0)
  if(length(m_gene_pos) > 0) {
    dat_gene <- tcs_iN_df[which(tcs_iN_df$Cell %in% colnames(dat_ref_tc_count[, m_gene_pos, drop = F])), ] %>%
      mutate(Gene = marker_genes_exo[i]) %>%
      mutate(TPM = as.numeric(dat_ref_tc_count[marker_genes_exo[i], which(colnames(dat_ref_tc_count) %in% .$Cell)])) %>%
      mutate(Gene = gsub("_EXO", "", .$Gene))
    final_dat_exo <- rbind(final_dat_exo, dat_gene)}}
final_dat_exo$TPM = log2(final_dat_exo$TPM + 1)

#plot
ggplot() +
  geom_point(data = tcs_iN_df, aes(UMAP1, UMAP2), size = 1.3, color = "grey40") +
  geom_point(data = final_dat_exo, aes(UMAP1, UMAP2, color = TPM), 
             size = 2) +
  scale_color_gradientn(colours = c("grey40", "white", "#ef9e45", "#f25c35", "#ee2c2c")) +
  theme_blank() +
  facet_wrap(~Gene, nrow = 4) +
  theme(strip.text = element_text(size = 20, color = "grey10"),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.8, colour = "grey10"),
        axis.ticks = element_line(size = 0.8, color = "grey10"),
        legend.position = "bottom")

####################################################################################################
############################################ FIGURE 4D #############################################
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

#create color palettes
my_pal_tc_monoc1 <- c("FIB" = "#033782",
                      "Ci1" = "#8b008b",
                      "Ci2" = "#c71585", 
                      "TFi1" = "#ee2c2c", 
                      "TFi2" = "#ee2c2c",
                      "TFi3" = "#ee7600", 
                      "TFi4" = "#ee7600", 
                      "a" = "white",
                      "b" = "white",
                      "c" = "#fbfcfe",
                      "d" = "#fbfcfe",
                      "e" = "#f6ebf9",
                      "f" = "#f8f2ea")

my_pal_tc_monoc2 <- c("FIB" = "#011736",
                      "Ci1" = "#400040", 
                      "Ci2" = "#7a0d52", 
                      "TFi1" = "#a11d1d", 
                      "TFi2" = "#a11d1d", 
                      "TFi3" = "#a15000", 
                      "TFi4" = "#a15000",
                      "a" = "#8b008b",
                      "b" = "#ec7600",
                      "c" = "#fbfcfe",
                      "d" = "#fbfcfe",
                      "e" = "#f6ebf9",
                      "f" = "#f8f2ea")

#create background
rect1 <- data.frame(xstart = c(-Inf, -5.5, 1, 1),
                    xend = c(-5.5, 1, +Inf, +Inf),
                    ystart = c(-Inf, -Inf, -1.1, -1.1),
                    yend = c(+Inf, +Inf, +Inf, -Inf),
                    color = c("c", "d", "e", "f"))

rect2 <- data.frame(xstart = c(1.1, 1.1),
                    xend = c(5.3, 5.3),
                    ystart = c(-1, -1.2),
                    yend = c(7.4, -9.2), 
                    color = c("a", "b"))

#plot with sample information
ggplot() +
  geom_rect(data = rect1, aes(xmin = xstart, 
                              xmax = xend, 
                              ymin = ystart, 
                              ymax = yend, 
                              color = color, 
                              fill = color), 
            size = 0.8) +
  geom_rect(data = rect2, aes(xmin = xstart, 
                              xmax = xend, 
                              ymin = ystart, 
                              ymax = yend, 
                              color = color, 
                              fill = color), 
            size = 0.8) +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"), 
               size = 1.2, 
               linetype = "solid", 
               na.rm = TRUE, 
               data = edge_df_iN) +
  geom_point(data = iNM_coord[sample(nrow(iNM_coord)), ], aes(V1, 
                                                              V2, 
                                                              fill = Sample, 
                                                              color = Sample), 
             size = 7, 
             shape = 21) +
  scale_fill_manual(values = my_pal_tc_monoc1) +
  scale_color_manual(values = my_pal_tc_monoc2) +
  theme_jo() +
  theme(panel.border = element_rect(fill = NA, 
                                    size = 1.5, 
                                    color = "grey10"))

#plot pseudotime information
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               size = 1.2, 
               linetype = "solid", 
               na.rm = TRUE, 
               data = edge_df_iN, 
               color = "grey20") +
  geom_point(data = iNM_coord, aes(V1, 
                                   V2, 
                                   fill = Pseudotime), 
             size = 12, 
             shape = 21) +
  scale_fill_gradient2(low = "#60b4b4", 
                       mid = "#004d81", 
                       high = "#f69c0d", 
                       midpoint = max(iNM_coord$Pseudotime/2)) +
  theme_jo() +
  theme(panel.border = element_rect(fill = NA, 
                                    size = 3, 
                                    color = "grey10"),
        axis.ticks = element_blank())

####################################################################################################
############################################ FIGURE 4E #############################################
####################################################################################################

genes <- readat(file.path(my_path, "Data", "gene_class.csv")) %>%
  filter(Class == "protein_coding") %>%
  .$Gene %>%
  .[which(. %in% expressed_genes)]

#BEAM
BEAM_res <- BEAM(iNM[genes, ], branch_point = 1, cores = 4)
BEAM_res_ord <- BEAM_res[order(BEAM_res$qval), ] %>%
  select(c("gene_short_name", "pval", "qval"))

#produce heatmap
beam_tc <- plot_genes_branched_heatmap(iNM[row.names(subset(BEAM_res_ord, qval < 1e-10)),], 
                                       norm_method = "log", 
                                       return_heatmap = T,
                                       num_clusters = 3)
beam_mat <- as.data.frame(beam_tc[[3]])

#plot heatmap
ph <- pheatmap(beam_mat, cluster_cols = F, cutree_rows = 4, clustering_method = "ward.D", fontsize_row = 1)

####################################################################################################
############################################ FIGURE 4F #############################################
####################################################################################################

#differential gene expression
dat_tc_ref <- readat(file.path(my_path, "Data", "tc_mat_tpm_qc_log2.csv"))
pt_dat <- pt_iN$data
my_states <- factor(pt_dat$State) 
names(my_states) <- pt_dat$sample_name
my_states <- my_states[names(tcs_iN@ident)]
tcs_iN@ident <- my_states

#find branch-specific marker genes
branch1_markers <- FindMarkers(tcs_iN, ident.1 = c(5,6,7,8,9), 
                               ident.2 = 4, 
                               logfc.threshold = 0.01, 
                               min.pct = 0.01, 
                               only.pos = T)
branch2_markers <- FindMarkers(tcs_iN,
                               ident.1 = 4, 
                               ident.2 = c(5,6,7,8,9), 
                               logfc.threshold = 0.01, 
                               min.pct = 0.01, 
                               only.pos = T)

#pvalue threshold of 0.05
branch1_markers_thres = branch1_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val <= 0.05)

branch2_markers_thres = branch2_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val <= 0.05)

#fold change threshold of 2
dat_split <- split(iNM_coord, iNM_coord$Branch)
dat_split_pos <- lapply(dat_split, function(X) rowMeans(subset(dat_tc_ref, select = X$Cell)))
dat_split_neg <- lapply(dat_split, function(X) rowMeans(dat_tc_ref[ ,-(which(names(dat_tc_ref) %in% X$Cell))]))
dat_fc <- as.data.frame(mapply("/", dat_split_pos, dat_split_neg, SIMPLIFY = T))
fold_branch1 <- dat_fc[2] %>%
  rownames_to_column(var = "gene") %>%
  filter(B > 2 & gene %in% branch1_markers_thres$gene)
fold_branch2 <- dat_fc[3] %>%
  rownames_to_column(var = "gene") %>%
  filter(C > 2 & gene %in% branch2_markers_thres$gene)

#final branch-specific marker tables
branch1_markers_final = branch1_markers_thres %>%
  filter(gene %in% fold_branch1$gene)
branch2_markers_final = branch2_markers_thres %>%
  filter(gene %in% fold_branch2$gene)

####################################################################################################
############################################ FIGURE 4G #############################################
####################################################################################################

#define exogenous TFs
marker_genes_exo <- c("ASCL1_EXO", "DLX1_EXO", "DLX2_EXO", "FEV_EXO", "FOXA2_EXO", "FOXP2_EXO", "ISL1_EXO", "LHX2_EXO",
                      "NEUROD1_EXO", "NEUROG2_EXO", "NR2F1_EXO", "NR2F2_EXO", "NR4A2_EXO", "OLIG2_EXO", "OTX2_EXO", 
                      "PAX6_EXO", "PITX3_EXO", "POU3F2_EXO", "TLX3_EXO", "ZIC1_EXO")

#create data frame
pt_dat <- pt_iN$data %>%
  mutate(State = as.numeric(as.character(State)),
         Branch = chartr("123456789", "AAACBBBBB", State))
pt_branch <- pt_dat %>%
  filter(Branch %in% c("B", "C")) %>%
  select(c("sample_name","Branch"))
pt_split <- pt_branch %>%
  split(.$Branch)
dat_split <- lapply(pt_split, function(X) dat_tc_ref[grep("_EXO", row.names(dat_tc_ref)), X$sample_name])

#perform fisher tests
my_fisher_pos <- matrix(nrow = 2, ncol = 0)

for (i in 1:length(marker_genes_exo)){
  exo <- marker_genes_exo[i]
  exo_pos <- length(which(dat_tc_ref[exo, pt_branch$sample_name] > 0))
  num_cells <- length(row.names(pt_branch))
  # Perform Fisher exact test
  x_list <- lapply(dat_split, function(x) length(which(x[exo, ] > 0)))
  m_list <- mapply(function(x, y) length(colnames(x))-y, dat_split, x_list, SIMPLIFY = F)
  n_list <- lapply(x_list, function(x) exo_pos - x)
  k_list <- mapply(function(x, y) num_cells - length(colnames(x)) - y, dat_split, n_list, SIMPLIFY = F)
  my_fisher <- mapply(function(a, b, c, d) fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
                     x_list, m_list, n_list, k_list, SIMPLIFY = F)
  my_p <- do.call("rbind", (lapply(my_fisher, function(x) x$p.value)))
  colnames(my_p) = exo
  my_fisher_pos = cbind(my_fisher_pos, my_p)
  print(i)
}

#modify fisher results
my_fisher_final <- as.data.frame(my_fisher_pos) %>%
  mutate(Branch = c("B", "C")) %>%
  melt(id = "Branch")
my_fisher_final$value <- -log10(my_fisher_final$value)

#create final data frame
final_branch <- iNM_coord %>%
  filter(Branch %in% c("B", "C")) %>%
  right_join(my_fisher_final, by = "Branch") %>%
  mutate(EXO = gsub("_EXO", "", variable))

#create color annotation
my_vec <- vector(mode = "character", length = length(row.names(final_branch)))
my_vec[which(final_branch$value >= 0 & final_branch$value <= 1.30103)] <- "A"
my_vec[which(final_branch$value > 1.30103 & final_branch$value <= 3)] <- "B"
my_vec[which(final_branch$value > 3 & final_branch$value <= 5)] <- "C"
my_vec[which(final_branch$value > 5)] <- "D"
final_branch$Color <- my_vec

#create color palette
my_pal_sig <- c("A" = "white", 
                "B" = "#f9d300", 
                "C" = "#fb7c00", 
                "D" = "#ed2323")

#plot
ggplot()+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"), 
               size = 1, 
               linetype = "solid", 
               na.rm = TRUE, 
               data = edge_df_iN) +
  geom_point(data = iNM_coord, aes(V1, V2), 
             size = 4.5, 
             shape = 21, 
             fill = "white", 
             color = "grey30") +
  geom_point(data = final_branch, aes(V1, V2, fill = Color), 
             color = "grey30", 
             shape = 21, 
             size = 4.5) +
  scale_fill_manual(values = my_pal_sig) +
  theme_blank() +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(size = 0.5)) +
  facet_wrap(~EXO, nrow = 4, scales = "free")

####################################################################################################
############################################ FIGURE 4H #############################################
####################################################################################################

#read IF Validation data
dat_IF_valid <- readex(file.path(my_path, "Data", "branch_valid.xlsx"), 3) %>%
  mutate(Sample = factor(Sample, levels = c("Pool", "Branch2", "Unsig")),
         Day = factor(Day, levels = c("9dpi", "21dpi")),
         Marker = factor(Marker, levels = c("TUBB3", "MAP2")))

#calculate mean and stdev
dat_IF_valid_sum <- summarise(group_by(dat_IF_valid, ID), Mean = mean(Count), STDEV = sd(Count)) %>%
  inner_join(dat_IF_valid, by = "ID") %>%
  select(-c(Exp, Count)) %>%
  mutate(Sample = factor(Sample, levels = c("Pool", "Branch2", "Unsig")),
         Day = factor(Day, levels = c("9dpi", "21dpi")),
         Marker = factor(Marker, levels = c("TUBB3", "MAP2"))) %>%
  unique %>%
  split(.$Marker)
dat_IF_valid <- split(dat_IF_valid, dat_IF_valid$Marker)

#create color palette
my_pal_IF_T1 <- c("Pool" = "grey40", "Branch2" = "#ed7600", "Unsig" = "#8a008a")
my_pal_IF_T2 <- c("Pool" = "grey20", "Branch2" = "#9e4e00", "Unsig" = "#3e003e")
my_pal_IF_T3 <- c("1" = 21, "2" = 21, "3" = 21, "4" = 21, "5" = 21)

#plot 9 dpi
ggplot(dat_IF_valid_sum[[1]], aes(Day, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean - STDEV, ymax = Mean + STDEV, color = Sample),
                stat = "identity", position = position_dodge(0.95), width = 0.7, size = 0.8) +
  geom_bar(stat = "identity", position = position_dodge(0.95), width = 0.8, size = 0.8) +
  geom_point(data = dat_IF_valid[[1]], aes(Day, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 4, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_IF_T1) +
  scale_color_manual(values = my_pal_IF_T2) +
  scale_shape_manual(values = my_pal_IF_T3) +
  theme_jo()

#plot 21 dpi
ggplot(dat_IF_valid_sum[[2]], aes(Day, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean - STDEV, ymax = Mean + STDEV, color = Sample),
                stat = "identity", position = position_dodge(0.95), width = 0.7, size = 0.8) +
  geom_bar(stat = "identity", position = position_dodge(0.95), width = 0.8, size = 0.8) +
  geom_point(data = dat_IF_valid[[2]], aes(Day, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 4, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_IF_T1) +
  scale_color_manual(values = my_pal_IF_T2) +
  scale_shape_manual(values = my_pal_IF_T3) +
  theme_jo()

####################################################################################################
############################################ FIGURE 4I #############################################
####################################################################################################

#read IF Validation data
dat_qpcr_valid <- readex(file.path(my_path, "Data", "branch_valid.xlsx"), 7) %>%
  mutate(Sample = factor(Sample, levels = c("Pool", "Branch2", "Unsig")),
         Day = factor(Day, levels = c("9dpi", "21dpi")),
         Marker = factor(Marker, levels = c("MAP2", "NCAM", "NEUN", "SYN1", "GLUT", "GABA", "TH", "CHAT", "VIM", "SNAI", "AREG", "PTHR1")))

#calculate mean and stdev
dat_qpcr_valid_sum <- summarise(group_by(dat_qpcr_valid, ID), Mean = mean(Count), STDEV = sd(Count)) %>%
  inner_join(dat_qpcr_valid, by = "ID") %>%
  select(-c(Exp, Count)) %>%
  mutate(Sample = factor(Sample, levels = c("Pool", "Branch2", "Unsig")),
         Day = factor(Day, levels = c("9dpi", "21dpi")),
         Marker = factor(Marker, levels = c("MAP2", "NCAM", "NEUN", "SYN1", "GLUT", "GABA", "TH", "CHAT", "VIM", "SNAI", "AREG", "PTHR1"))) %>%
  unique 

#pan markers
dat_pan <- dat_qpcr_valid %>%
  filter(Marker %in% c("MAP2", "NCAM", "NEUN", "SYN1")) %>%
  split(.$Day)

dat_pan_sum <- dat_qpcr_valid_sum %>%
  filter(Marker %in% c("MAP2", "NCAM", "NEUN", "SYN1")) %>%
  split(.$Day)

#create color palette
my_pal_qpcr_T1 <- c("Pool" = "grey40", "Branch2" = "#ed7600", "Unsig" = "#8a008a")
my_pal_qpcr_T2 <- c("Pool" = "grey20", "Branch2" = "#9e4e00", "Unsig" = "#3e003e")
my_pal_qpcr_T3 <- c("1" = 21, "2" = 21, "3" = 21, "4" = 21)

#plot 9dpi
ggplot(dat_pan_sum[[1]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_pan[[1]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 16)) +
  theme_jo()

#plot 21dpi
ggplot(dat_pan_sum[[2]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_pan[[2]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 16)) +
  theme_jo()

####################################################################################################
############################################ FIGURE 4J #############################################
####################################################################################################

#subtype markers
dat_sub <- dat_qpcr_valid %>%
  filter(Marker %in% c("GLUT", "GABA", "TH", "CHAT")) %>%
  split(.$Day)

dat_sub_sum <- dat_qpcr_valid_sum %>%
  filter(Marker %in% c("GLUT", "GABA", "TH", "CHAT")) %>%
  split(.$Day)

#plot 9dpi
ggplot(dat_sub_sum[[1]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_sub[[1]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 20)) +
  theme_jo()

#plot 21dpi
ggplot(dat_sub_sum[[2]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_sub[[2]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 20)) +
  theme_jo()

####################################################################################################
############################################ FIGURE 4K #############################################
####################################################################################################

#fibroblast markers
dat_fib <- dat_qpcr_valid %>%
  filter(Marker %in%  c("VIM", "SNAI")) %>%
  split(.$Day)
dat_fib_sum <- dat_qpcr_valid_sum %>%
  filter(Marker %in%  c("VIM", "SNAI")) %>%
  split(.$Day)

#plot 9dpi
ggplot(dat_fib_sum[[1]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_fib[[1]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2)) +
  theme_jo()

#plot 21dpi
ggplot(dat_fib_sum[[2]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_fib[[2]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.2)) +
  theme_jo()

####################################################################################################
############################################ FIGURE 4K #############################################
####################################################################################################

#vascular markers
dat_vasc <- dat_qpcr_valid %>%
  filter(Marker %in% c("AREG", "PTHR1")) %>%
  split(.$Day)
dat_vasc_sum <- dat_qpcr_valid_sum %>%
  filter(Marker %in% c("AREG", "PTHR1")) %>%
  split(.$Day)

#plot 9dpi
ggplot(dat_vasc_sum[[1]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_vasc[[1]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 18)) +
  theme_jo()  

#plot 9dpi
ggplot(dat_vasc_sum[[2]], aes(Marker, Mean, fill = Sample, color = Sample)) +
  geom_errorbar(aes(ymin = Mean-STDEV, ymax = Mean+STDEV, color = Sample), 
                position = position_dodge(0.9), width = 0.5, size = 0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75, size = 0.6) +
  geom_point(data = dat_vasc[[2]], aes(Marker, Count, fill = Sample, shape = as.factor(Exp)), 
             color = "grey10", size = 5, position = position_dodge(0.95), stroke = 0.6) +
  scale_fill_manual(values = my_pal_qpcr_T1) +
  scale_color_manual(values = my_pal_qpcr_T2) +
  scale_shape_manual(values = my_pal_qpcr_T3) +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 18)) +
  theme_jo()