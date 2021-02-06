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
########################################### FIGURE S3A #############################################
####################################################################################################

#read Seurat object
iNS_10x <- get(load(file.path(my_path, "Data", "Seurat_iN_10x.Rmd")))

#create data frame for ggplot
iNS_10x_df <- as.data.frame(iNS_10x@dr$umap@cell.embeddings) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(Sample = factor(unlist(lapply(strsplit(Cell, "[.]"), function(X) X[[1]])), levels = c("Ci", "TFi")),
         Cluster = factor(iNS_10x@meta.data$res.0.7, levels = c(1,2,3,4,5,6,7,8,9)))

#read 10x exo data
dat_exo <- readat(file.path(my_path, "Data", "dat_reference_10x_exo.csv"))
dat_exo <- log2(dat_exo + 1)

#define exogenous and endogenous genes
marker_genes_exo <- c("ASCL1_EXO","DLX1_EXO","DLX2_EXO","FEV_EXO","FOXA2_EXO",
                     "ISL1_EXO","NEUROD1_EXO","NR2F1_EXO","NR2F2_EXO","NR4A2_EXO",
                     "OLIG2_EXO","PAX6_EXO","PITX3_EXO","POU3F2_EXO","ZIC1_EXO")

marker_genes_endo <- gsub("_EXO", "", marker_genes_exo)

#extract TPM of EXO marker genes
final_dat_exo <- matrix(ncol = 7, nrow = 0)
for (i in 1:length(marker_genes_exo)) {
  m_gene_pos <- which(dat_exo[marker_genes_exo[i], ] > 0)
  if(length(m_gene_pos) > 0) {
    dat_gene <- iNS_10x_df[which(iNS_10x_df$Cell %in% colnames(dat_exo[, m_gene_pos, drop = F])), ] %>%
      mutate(Gene = marker_genes_exo[i]) %>%
      mutate(TPM = as.numeric(dat_exo[marker_genes_exo[i], which(colnames(dat_exo) %in% .$Cell)])) %>%
      mutate(Gene = gsub("_EXO", "", .$Gene))
    final_dat_exo <- rbind(final_dat_exo, dat_gene)}}

# Extract TPM of ENDO marker genes
final_dat_endo <- matrix(ncol = 7, nrow = 0)

for (i in 1:length(marker_genes_endo)) {
  m_gene_pos <- which(dat_exo[marker_genes_endo[i], ] > 0)
  if(length(m_gene_pos) > 0) {
    dat_gene <- iNS_10x_df[which(iNS_10x_df$Cell %in% colnames(dat_exo[, m_gene_pos, drop = F])), ] %>%
      mutate(Gene = marker_genes_endo[i]) %>%
      mutate(TPM = as.numeric(dat_exo[marker_genes_endo[i], which(colnames(dat_exo) %in% .$Cell)]))
    final_dat_endo = rbind(final_dat_endo, dat_gene)}}

#Plot EXO data in facet grid
ggplot() +
  geom_point(data = iNS_10x_df, aes(UMAP1, UMAP2), size = 1.3, color = "grey60", alpha = 0.6) +
  geom_point(data = final_dat_exo, aes(UMAP1, UMAP2, color = TPM), 
             size = 1.8, alpha = 1) +
  scale_color_gradient2(limits = c(min(final_dat_exo$TPM), max(final_dat_exo$TPM)), 
                        low = "#9d9d9d", 
                        mid = "white", 
                        high = "#ac0002",
                        midpoint=(max(final_dat_exo$TPM)/2)) +
  theme_blank() +
  facet_wrap(~Gene, ncol = 5, nrow = 3) +
  theme(strip.text = element_text(size = 20, color = "grey10"),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.8, colour = "grey10"),
        axis.ticks = element_line(size = 0.8, color = "grey10"),
        legend.position = "bottom")

#Plot ENDO data in facet grid
ggplot() +
  geom_point(data = iNS_10x_df, aes(UMAP1, UMAP2), size = 1.3, color = "grey60", alpha = 0.6) +
  geom_point(data = final_dat_endo, aes(UMAP1, UMAP2, color = TPM), 
             size = 1.8) +
  scale_color_gradient2(limits = c(0, max(final_dat_endo$TPM)), 
                        low = "#9d9d9d", 
                        mid = "white", 
                        high = "#1b75ba",
                        midpoint = (max(final_dat_endo$TPM)/2)) +
  theme_blank()+
  facet_wrap(~Gene, ncol = 5, nrow = 3) +
  theme(strip.text = element_text(size = 20, color = "grey10"),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.8, colour = "grey10"),
        axis.ticks = element_line(size = 0.8, color = "grey10"),
        legend.position = "bottom")

####################################################################################################
########################################### FIGURE S3B #############################################
####################################################################################################

#read data
read_list <- list()
for (i in 1:7) {
  read_list[[i]] <- readex(file.path(my_path, "Data", "bowtie_genome_browser.xlsx"), i)}
incept <- read_list[[5]]

#define end of genome browser
gp_end <- max(read_list[[7]]$V8+100)

#create color palette
my_pal_genome <- c("ASCL1_EXO" = "#df5e00", "DLX1_EXO" = "#d0352a", "DLX2_EXO" = "#bf8c00", "FEV_EXO" = "#a3524d", 
                  "FOXA2_EXO" = "#710aff", "FOXP2_EXO" = "#9f95ff", "ISL1_EXO" = "#326ddd", "LHX2_EXO" = "#7fffce", 
                  "NEUROD1_EXO" = "#5da8b5", "NGN2_EXO" = "#24ab2b", "NR2F1_EXO" = "#ff9c00", "NR2F2_EXO" = "#bdba00",
                  "NR4A2_EXO" = "#ff9a9c", "OLIG2_EXO" = "#a936af", "OTX2_EXO" = "#ce746f", "PAX6_EXO" = "#ffd600", 
                  "PITX3_EXO" = "#dfc496", "POU3F2_EXO" = "#db3082", "TLX3_EXO" = "#8883c9", "ZIC1_EXO" = "#906a91",
                  "3" = "#de8800", "4" = "#5c5c5c")

#POOL
dat_read_POOL <- read_list[[7]] %>%
  split(.$V1) %>%
  lapply(function(y) data.frame(Count = colSums(genome.browser(y, read_list[[7]])), 
                                Position = c(1:gp_end))) %>%
  do.call("cbind",.) %>%
  select(-ends_with("Position")) %>%
  mutate(Position = c(1:gp_end)) %>%
  melt(id = "Position") %>%
  mutate(variable = gsub(".Count", "", variable), 
         intercept1 = rep(incept$Intercept1, each = gp_end),
         intercept2 = rep(incept$Intercept2, each = gp_end)) %>%
  split(.$variable) %>%
  lapply(function(y) {y$Color.i1 = ifelse(y[unique(y$intercept1), 3] > 0, 3, 4); y}) %>%
  lapply(function(y) {y$Color.i2 = ifelse(y[unique(y$intercept2), 3] > 0, 3, 4); y}) %>%
  do.call("rbind",.) %>%
  mutate(arrow1 = max(value),
         arrow2 = max(value)/2)

#plot
ggplot(dat_read_POOL, aes(Position, value)) +
  geom_density(aes(fill = as.factor(variable), color = as.factor(variable)), stat = "identity") +
  geom_segment(aes(intercept1, arrow1, xend = intercept1, yend = arrow2, color = as.factor(Color.i1)), size = 0.4) +
  geom_segment(aes(intercept2, arrow1, xend = intercept2, yend = arrow2, color = as.factor(Color.i2)), size = 0.4) +
  scale_color_manual(values = my_pal_genome) +
  scale_fill_manual(values = my_pal_genome) +
  xlim(c(900, gp_end)) +
  theme_blank() +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(size = 0.7)) +
  facet_wrap(~variable, ncol = 1) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 0.5, color = "black") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 0.5, color = "black")

####################################################################################################
########################################### FIGURE S3C #############################################
####################################################################################################

#set the directory to the main folder and read all files in all subfolders
read_list <- list()
for (i in 1:4) {
  read_list[[i]] <- readex(file.path(my_path, "Data", "kallisto_convertseq.xlsx"), i)}
sf <- list.files(file.path(my_path, "Data", "g24_kallisto_full_backloading"), 
                recursive = T, 
                pattern = ("abundance.tsv"), 
                full = T)
all_files <- pbapply::pblapply(sf, fread)

#merge files
all_tpm <- lapply(all_files, function(y) as.data.frame(y$tpm)) %>%
  do.call("cbind", .)
names(all_tpm) <- read_list[[2]]$NAMES
s1 <- strsplit(as.character(all_files[[1]]$target_id), "[.]")
s2 <- sapply(1:length(s1), function(X) s1[[X]][1])
all_tpm$s2 <- s2
m1 <- merge(all_tpm, read_list[[1]], by.x = "s2", by.y = "ID")
m2 <- m1[,-1]

#summarize 
dat_tpm <- plyr::ddply(m2, "Gene", plyr::numcolwise(sum), .progress = "text") %>%
  filter(Gene %in% c(read_list[[4]]$Genes_EXO, read_list[[4]]$Genes)) %>%
  column_to_rownames(var = "Gene") %>%
  as.data.frame %>%
  .[c(read_list[[4]]$Genes_EXO, read_list[[4]]$Genes), read_list[[3]]$LABEL] %>%
  mutate_all(function(x) log2(x+1)) %>%
  as.matrix
row.names(dat_tpm) <- rep(read_list[[4]]$Genes, times = 2)

#split exogenous and endogenous
dat_exo <- dat_tpm[1:20, ]
dat_endo <- dat_tpm[21:40, ]

#plot heatmaps
pheatmap(dat_exo, cluster_rows = FALSE, cluster_cols = F, border_color = "grey40", fontsize = 15, 
         show_colnames = F, color = colorRampPalette(c("#4d4d4d", "white", "#ad0002"))(100), legend = T)

pheatmap(dat_endo, cluster_rows = FALSE, cluster_cols = F, border_color = "grey40", fontsize = 15, 
         show_colnames = F, color = colorRampPalette(c("#4d4d4d", "white", "#0d6dba"))(100), legend = T)