#check if pacman package is loaded and install if not
if (!require(pacman)) {
  install.packages("pacman")
  if(!require(pacman)) stop("Package pacman not found")
}

#set data path
my_path <- "D:/Projects/Single_Cell_Convert_Seq"

#load libraries and custom functions
pacman::p_load(rio, tidyverse, readxl, reshape, edgeR)
source(file.path(my_path, "R", "utilities.R"))

####################################################################################################
############################################ FIGURE S1A ############################################
####################################################################################################

#read data
read_list <- list()
for (i in 1:4) {
  read_list[[i]] <- readex(file.path(my_path, "Data", "ips_to_neuron.xlsx"), i)}
dat_ips <- read_list[[1]] %>%
  column_to_rownames(var = "Gene") %>%
  as.data.frame()
dat_nervous <- read_list[[2]]
dat_cand <- read_list[[3]] %>%
  mutate(p.Gene = paste0("p1@", .$Gene))
dat_tf <- read_list[[4]]

#define groups and subset data
group <- dat_ips %>% 
  colnames %>% 
  stringr::str_split("_") %>% 
  sapply(., function(X) X[3]) %>% 
  .[1:12]
f_data <- dat_ips[rowSums(dat_ips) > 5, 1:12]

#combine data and group; Calculate dispersions; Fit dispersion
my_data <- f_data %>% DGEList(counts = ., group = group) %>% estimateCommonDisp %>% estimateTagwiseDisp
design <- model.matrix(~group)
colnames(design) = unique(group)
fit <- glmFit(my_data, design)

#ANOVA-like analysis
lrt <- glmLRT(fit, coef = 2:4)
lrt$table$FDR <- p.adjust(lrt$table$PValue, method = "BH")
de_table <- lrt$table

#set a cutoff, add day0
sig_data <- de_table %>%
  rownames_to_column %>%
  filter(FDR <= 0.05) %>% 
  select(1:4) %>% 
  setNames(gsub("logFC.", "", names(.))) %>% 
  mutate(day00 = 0) %>% 
  column_to_rownames(var = "rowname") %>%
  select("day00", "day06", "day12", "day18") %>%
  t %>% 
  as.data.frame

#select TFs
tfs <- dat_tf %>% 
  mutate(p.Gene = paste0("p1@", .$Gene))
sig_data_tf <- sig_data[, which(colnames(sig_data) %in% tfs$p.Gene)] %>%
  t %>% 
  melt %>% 
  full_join(dat_cand, by = c("X1" = "p.Gene")) %>% 
  select(c("X1", "X2", "value", "Color"))
sig_data_tf = sig_data_tf[-(which(is.na(sig_data_tf$value) == "TRUE")), ]
sig_data_tf[is.na(sig_data_tf)] = "Z"

#create color palette
my_pal_ips <- c("Z" = "grey90", "A" = "#1b75bb", "B" = "#f6b520", "C" = "#bf282c", "D" = "#545453")

#ploy data
ggplot(sig_data_tf, aes(X2, value, group = X1, color = as.factor(Color), fill = as.factor(Color))) +
  geom_line(alpha = 0.5, lwd = 1, lineend = "round") +
  geom_point(data = sig_data_tf[-(which(sig_data_tf$Color == "Z")), ], size = 4, pch = 19) +
  geom_point(data = sig_data_tf[-(which(sig_data_tf$Color == "Z")), ], size = 2, pch = 19, color = "white") +
  geom_line(data = sig_data_tf[-(which(sig_data_tf$Color == "Z")), ], lwd = 1, lineend = "round") +
  scale_color_manual(values = my_pal_ips) +
  scale_fill_manual(values = my_pal_ips) +
  scale_x_discrete(expand = c(0.1, 0.1), labels = c("Day0", "Day6", "Day12", "Day18")) +
  labs(x = "Time-point", y = "Log2(FC vs. Day0)") +
  geom_hline(aes(yintercept = 0), color = "grey20", lty = 2, size = 1) +
  theme_jo() +
  theme(legend.position = "none")

####################################################################################################
############################################ FIGURE S1C ############################################
####################################################################################################
#run model
my_per <- c(3, 5, 7)
final_frame <- matrix(nrow = 0, ncol = 2)
for (j in 1:length(my_per)) {
  #create list
  my_list = list()  
  for(l in seq_len(1000)) {
    mat <- matrix(data = 0, nrow = 20, ncol = 1000)
    colnames(mat) = paste0("Cell", 1:1000)
    row.names(mat) = paste0("TF", 1:20)
    for (i in seq_len(20)){
      mat[i, ][sample(1000, 1000/my_per[j])] = 1
    }
    x <- colSums(mat)
    my_list[[l]] <- as.data.frame(table(x))
  }
  
  #create data frame for plotting
  my_frame <- my_list %>% 
    do.call("rbind", .) %>% 
    split(.$x) %>% 
    lapply(function(X) mean(X$Freq)) %>% 
    unlist %>% 
    as.data.frame %>%
    plyr::rename(replace = c("." = "Freq")) %>%
    rownames_to_column(var = "TF")
  my_freq <- freq(my_frame) %>%
    mutate(Group = LETTERS[j])
  
  final_frame <- rbind(final_frame, my_freq)
}

#color palette
my_pal <- c("A" = "#535353", "B" = "#C0292C", "C" = "grey70")

#plot data
ggplot(final_frame, aes(V1, group = Group, fill = Group, color = Group)) +
  geom_density(bw = 1, alpha = 0.4) +
  #scale_x_discrete(breaks = c(0, 5, 10, 15, 20)) +
  scale_color_manual(values = my_pal) +
  scale_fill_manual(values = my_pal) +
  theme_jo()