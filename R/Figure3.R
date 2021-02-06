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

####################################################################################################
############################################ FIGURE 3B #############################################
####################################################################################################

#read data
read_list <- list()
for (i in 1:6) {
  read_list[[i]] = readex(file.path(my_path, "Data", "bowtie_genome_browser.xlsx"), i)}

#filter for selected TFs
selected_TFs <- c("DLX1_EXO", "FOXA2_EXO", "ISL1_EXO", "NEUROD1_EXO", "NGN2_EXO", 
                 "NR2F1_EXO", "OLIG2_EXO", "PAX6_EXO", "POU3F2_EXO", "ZIC1_EXO")
read_list <- lapply(read_list, function(X) filter(X, V1 %in% selected_TFs))
incept <- read_list[[5]]

#define end of genome browser
gp_end <- max(read_list[[6]]$V3 + 100)

#DLX1
dat_read_DLX1 <- read_list[[6]] %>%
  filter(ID %in% c("A12", "C10", "E08", "G06")) %>%
  select(-ID) %>%
  split(.$V1) %>%
  lapply(function(y) data.frame(Count = colSums(genome.browser(y, read_list[[6]])), 
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

#COMB1
dat_read_COMB1 <- read_list[[6]] %>%
  filter(ID %in% c("H06", "H07")) %>%
  select(-ID) %>%
  split(.$V1) %>%
  lapply(function(y) data.frame(Count = colSums(genome.browser(y, read_list[[6]])), 
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

#COMB2
dat_read_COMB2 <- read_list[[6]] %>%
  filter(ID %in% c("H08", "H09")) %>%
  select(-ID) %>%
  split(.$V1) %>%
  lapply(function(y) data.frame(Count = colSums(genome.browser(y, read_list[[6]])), 
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

#Ci
dat_read_Ci = as.data.frame(read_list[[4]]) %>%
  mutate(intercept1 = rep(incept$Intercept1, each = 4035),
         intercept2 = rep(incept$Intercept2, each = 4035)) %>%
  split(.$V1) %>%
  lapply(function(y) {y$Color.i1 = ifelse(y[unique(y$intercept1), 3] > 0, 3, 4); y}) %>%
  lapply(function(y) {y$Color.i2 = ifelse(y[unique(y$intercept2), 3] > 0, 3, 4); y}) %>%
  do.call("rbind",.) %>%
  mutate(Color = 2, arrow1 = 10, arrow2 = 5)

#create color palette
my_pal_genome = c("DLX1_EXO" = "#d0352a", 
                  "FOXA2_EXO" = "#949494", 
                  "ISL1_EXO" = "#949494", 
                  "NEUROD1_EXO" = "#949494", 
                  "NGN2_EXO" = "#949494",
                  "NR2F1_EXO" = "#949494", 
                  "OLIG2_EXO" = "#949494", 
                  "PAX6_EXO" = "#949494", 
                  "POU3F2_EXO" = "#949494", 
                  "ZIC1_EXO" = "#949494",
                  "3" = "#de8800", 
                  "4" = "#5c5c5c")

#plot DLX1
ggplot(dat_read_DLX1, aes(Position, value)) +
  geom_density(aes(fill = as.factor(variable), color = as.factor(variable)), stat = "identity") +
  geom_segment(aes(intercept1, arrow1, xend = intercept1, yend = arrow2, color = as.factor(Color.i1)), size = 0.4) +
  geom_segment(aes(intercept2, arrow1, xend = intercept2, yend = arrow2, color = as.factor(Color.i2)), size = 0.4) +
  scale_color_manual(values = my_pal_genome) +
  scale_fill_manual(values = my_pal_genome) +
  xlim(c(900, gp_end)) +
  theme_blank() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~variable, ncol = 1) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 0.5, color = "black") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 0.5, color = "black")

#create color palette
my_pal_genome = c("DLX1_EXO" = "#d0352a", 
                  "FOXA2_EXO" = "#a3524d", 
                  "ISL1_EXO" = "#326ddd", 
                  "NEUROD1_EXO" = "#5da8b5", 
                  "NGN2_EXO" = "#24ab2b",
                  "NR2F1_EXO" = "#949494", 
                  "OLIG2_EXO" = "#949494", 
                  "PAX6_EXO" = "#949494", 
                  "POU3F2_EXO" = "#949494", 
                  "ZIC1_EXO" = "#949494",
                  "3" = "#de8800", 
                  "4" = "#5c5c5c")
#plot COMB1
ggplot(dat_read_COMB1, aes(Position, value)) +
  geom_density(aes(fill = as.factor(variable), color = as.factor(variable)), stat = "identity") +
  geom_segment(aes(intercept1, arrow1, xend = intercept1, yend = arrow2, color = as.factor(Color.i1)), size = 0.4) +
  geom_segment(aes(intercept2, arrow1, xend = intercept2, yend = arrow2, color = as.factor(Color.i2)), size = 0.4) +
  scale_color_manual(values = my_pal_genome) +
  scale_fill_manual(values = my_pal_genome) +
  xlim(c(900, gp_end)) +
  theme_blank() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~variable, ncol = 1) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 0.3, color = "black") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 0.3, color = "black")

#create color palette
my_pal_genome = c("DLX1_EXO" = "#949494", 
                  "FOXA2_EXO" = "#949494", 
                  "ISL1_EXO" = "#949494", 
                  "NEUROD1_EXO" = "#949494", 
                  "NGN2_EXO" = "#949494",
                  "NR2F1_EXO" = "#ff9c00", 
                  "OLIG2_EXO" = "#a936af", 
                  "PAX6_EXO" = "#ffd600", 
                  "POU3F2_EXO" = "#db3082", 
                  "ZIC1_EXO" = "#906a91",
                  "3" = "#de8800", 
                  "4" = "#5c5c5c")
#plot COMB2
ggplot(dat_read_COMB2, aes(Position, value)) +
  geom_density(aes(fill = as.factor(variable), color = as.factor(variable)), stat = "identity") +
  geom_segment(aes(intercept1, arrow1, xend = intercept1, yend = arrow2, color = as.factor(Color.i1)), size = 0.4) +
  geom_segment(aes(intercept2, arrow1, xend = intercept2, yend = arrow2, color = as.factor(Color.i2)), size = 0.4) +
  scale_color_manual(values = my_pal_genome) +
  scale_fill_manual(values = my_pal_genome) +
  xlim(c(900, gp_end)) +
  theme_blank() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~variable, ncol = 1) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 0.3, color = "black") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 0.3, color = "black")

#plot Ci
ggplot(dat_read_Ci, aes(Position, value)) +
  geom_density(fill = "#5c5c5c", color = "#5c5c5c", stat = "identity") +
  geom_segment(aes(intercept1, arrow1, xend = intercept1, yend = arrow2, color = as.factor(Color.i1)), size = 0.4) +
  geom_segment(aes(intercept2, arrow1, xend = intercept2, yend = arrow2, color = as.factor(Color.i2)), size = 0.4) +
  scale_color_manual(values = my_pal_genome) +
  scale_fill_manual(values = my_pal_genome) +
  ylim(c(0, 10)) + 
  xlim(c(900, gp_end)) +
  theme_blank() +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~V1, ncol = 1) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 0.5, color = "black") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 0.5, color = "black")

####################################################################################################
############################################ FIGURE 3C #############################################
####################################################################################################

#read data
read_list <- list()
for (i in 1:2) {
  read_list[[i]] = readex(file.path(my_path, "Data", "bowtie_convertseq.xlsx"), i)}

dat_bowtie <- read_list[[1]] %>%
  column_to_rownames(var = "ID") %>%
  mutate_all(function(x) log2(x+1))

#plot
pheatmap(as.matrix(dat_bowtie), cluster_rows = F, cluster_cols = F, border_color = "grey40", fontsize = 15, show_colnames = F, 
         color = colorRampPalette(c("#4d4d4d", "white", "#ad0002"))(100), annotation_legend = F, legend = F)

####################################################################################################
############################################ FIGURE 3D #############################################
####################################################################################################

#set the directory to the main folder and read all files in all subfolders
#read data
read_list <- list()
for (i in 1:4) {
  read_list[[i]] <- readex(file.path(my_path, "Data", "kallisto_convertseq.xlsx"), i)}
sf <- list.files(file.path(my_path, "Data", "g24_kallisto_backloading"), recursive = T, pattern = ("abundance.tsv"), full = T)
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

#plot heatmap
pheatmap(dat_tpm, cluster_rows = FALSE, cluster_cols = F, border_color = "grey40", fontsize = 15, 
         show_colnames = F, color = colorRampPalette(c("#4d4d4d", "white", "#B72026"))(100), legend = T)

####################################################################################################
############################################ FIGURE 3E #############################################
####################################################################################################

#read data
dat_exendo <- readat(file.path(my_path, "Data", "backloading_exo_endo.csv"))
dat_exo <- dat_exendo %>%
  filter(Origin == "EXO")
dat_out <- dat_exendo %>%
  filter(Out == 4)

#create color palettes
my_pal_back1 <- c("#08406e", "#5e0002")
my_pal_back2 <- c("#0d6ebb", "#ac0002")

ggplot(dat_exendo) +
  geom_boxplot(aes(Gene, TPM, color = as.factor(Fill), fill = as.factor(Fill)), 
               outlier.shape = 15, outlier.size = 0.5, outlier.alpha = 0.7) +
  geom_point(data = dat_out, aes(Gene, TPM), fill = "#dc8101", color= "#8f5401", 
             shape = 21, size = 2, alpha = 1) + 
  stat_summary(data = dat_out, geom = "errorbar", fun.y = mean, 
               aes(Gene, TPM, group = Group, ymin = ..y.., ymax = ..y..), 
               color = "#dc8101", size = 0.7, alpha = 1) +
  scale_fill_manual(values = my_pal_back2) +
  scale_color_manual(values = my_pal_back1) +
  ylim(c(0, 8)) +
  theme_jo() +
  theme(panel.grid.major.y = element_line(color = "grey80"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 8))