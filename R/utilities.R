readat <- function(id) {
  dat = import(id)
  dat = dat %>% remove_rownames %>% column_to_rownames(var = names(dat)[1])
}

readex <- function(id, sheet) {
  dat = read_excel(id, sheet = sheet)
}

freq <- function(frame){
  mat = matrix(ncol = 1,nrow = 0)
  for (i in 1:length(row.names(frame))){
    x = data.frame(V1 = rep(frame[i, 1], times = frame[i, 2]))
    mat = rbind(mat, x)}
  return(mat)
}

genome.browser <- function(frame, dat){
  mat = matrix(0, length(row.names(frame)), max(dat$V3)+100)
  for (i in 1:length(row.names(frame))) {
    mat[i, c(frame$V2[i]:frame$V3[i])] = 1}
  return(mat)
}

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

theme_jo <- function(base_size = 12){
  theme_classic(base_size = base_size) +
    theme(
      legend.key = element_blank(), 
      strip.background = element_blank(),
      strip.text = element_text(size = 20, color = "grey10"),
      plot.title = element_blank(), 
      axis.title = element_text(size = 15, color = "grey10"), 
      axis.text = element_text(size = 15, color = "grey10"), 
      legend.title = element_blank(), 
      legend.text = element_blank(),
      axis.line = element_line(size = 0.8, color = "grey10"),
      axis.ticks = element_line(size = 0.8, color = "grey10"),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none")
}

theme_blank <- function(base_size = 12){
  theme_jo(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      strip.text = element_blank()
    )
}

darken <- function(color, factor = 1.5){
  col = col2rgb(color)
  col = col/factor
  col = rgb(t(col), maxColorValue = 255)
  names(col) = names(color)
  col
}

Feature.Jo <- function(Seurat, 
                       Genes, 
                       Overlay = FALSE,
                       Reduction = "UMAP",
                       Col1 = "#ffa017",
                       Col2 = "#df1900") {
  
  #darken function
  darken <- function(color, factor = 1.5){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue = 255)
    names(col) <- names(color)
    col
  }
  
  #create theme
  theme_feature <- function(base_size = 12){
    theme_blank(base_size = base_size) +
      theme(
        plot.title = element_text(size = 60),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 25),
        legend.position = "bottom",
        legend.text = element_text(),
        legend.title = element_text())
  }
  
  #create data frames with reduction coordinates
  if (Reduction == "UMAP") {
    
    df1 <- as.data.frame(Embeddings(Seurat[["umap"]]))
    
  } else if (Reduction == "TSNE") {
      
    df1 <- as.data.frame(Embeddings(Seurat[["tsne"]]))
    
    }
  
  
  Seurat.df <- df1 %>%
    rownames_to_column(var = "Cell")
  
  #check input Genes
  m <- setdiff(Genes, rownames(Seurat@assays$RNA))
  
  if (length(m) > 0) {
    stop(paste0("The genes ", m, " do not exist in Seurat@data"))
    
  } else {
    #log transform count data
    dat.ref <- Seurat@assays$RNA@data[which(rownames(Seurat@assays$RNA@data) %in% Genes), , drop = F]
    dat.ref <- log2(dat.ref + 1)
  }
  
  #create data frame with gene expression
  final.markers <- matrix(ncol = 6, nrow = 0)
  
  for (i in Genes) {
    
    if (length(which(dat.ref[i, ]>0))>0) {
      
      gene.pos <- which(dat.ref[i, ]>0)
      dat.gene <- Seurat.df %>% 
        filter(Cell %in% colnames(dat.ref[, gene.pos, drop = F])) %>%
        mutate(Gene = i,
               TPM = as.numeric(dat.ref[i, which(colnames(dat.ref) %in% .$Cell)]))
      final.markers <- rbind(final.markers, dat.gene)
    } else {
      print(paste0("Gene ", i, " is not expressed in any cell"))
    }
  }
  
  if (Overlay) {
    
    if (length(Genes) != 2) {
      stop("Overlay only works for 2 Genes")
    }
    
    exclude <- which(Seurat.df$Cell %in% final.markers$Cell)
    col.pal <- c(Col1, Col2)
    fill.pal <- darken(col.pal)
    
    if (Reduction == "UMAP") {
      
      print(ggplot() +
              geom_point(data = Seurat.df[-exclude, ], aes(UMAP1, UMAP2),
                         size = 5,
                         color = "grey70",
                         alpha = 0.5) +
              geom_point(data = final.markers[sample(rownames(final.markers)),], aes(UMAP1, UMAP2,
                                                                                     size = TPM,
                                                                                     color = Gene,
                                                                                     fill = Gene),
                         alpha = 0.7) +
              scale_color_manual(values = col.pal) +
              scale_fill_manual(values = fill.pal) +
              guides(size = guide_legend(title = "Log2(CPM)"),
                     color = guide_legend(title = "Genes"),
                     fill = FALSE) +
              scale_size(range = c(0, 7)) +
              theme_feature())
      
    } else if (Reduction == "TSNE") {
      
      print(ggplot() +
              geom_point(data = Seurat.df[-exclude, ], aes(tSNE_1, tSNE_2),
                         size = 5,
                         color = "grey70",
                         alpha = 0.5) +
              geom_point(data = final.markers[sample(rownames(final.markers)),], aes(tSNE_1, tSNE_2,
                                                                                     size = TPM,
                                                                                     color = Gene,
                                                                                     fill = Gene),
                         alpha = 0.7) +
              scale_color_manual(values = col.pal) +
              scale_fill_manual(values = fill.pal) +
              guides(size = guide_legend(title = "Log2(CPM)"),
                     color = guide_legend(title = "Genes"),
                     fill = FALSE) +
              scale_size(range = c(0, 7)) +
              theme_feature())
    } else {
      stop("Unknown reduction type")
    }
    
  } else {
    
    final.markers <- final.markers %>%
      mutate(Gene <- factor(Gene, levels = Genes)) %>%
      as.data.table %>%
      .[order(Gene, TPM)] 
    
    if (Reduction == "UMAP") {
      
      print(ggplot() +
              geom_point(data = Seurat.df, aes(UMAP_1, UMAP_2),
                         color = "#cdcdcb",
                         size = 1) +
              geom_point(data = final.markers, aes(UMAP_1, UMAP_2, color = TPM),
                         size = 1) +
              scale_color_gradientn("Log2(counts)", 
                                    colors = colorRampPalette(c("#cdcdcb", 
                                                                "#e7f5b6", 
                                                                "#3db2c4", 
                                                                "#0d3286",
                                                                "#043b5c"))(100)) +
              theme_feature() +
              facet_wrap(~Gene))
      
    } else if (Reduction == "TSNE") {
      
      print(ggplot() +
              geom_point(data = Seurat.df, aes(tSNE_1, tSNE_2),
                         color = "#cdcdcb",
                         size = 3) +
              geom_point(data = final.markers, aes(tSNE_1, tSNE_2, color = TPM),
                         size = 3) +
              scale_color_gradientn("Log2(counts)",
                                    colors = colorRampPalette(c("#cdcdcb", 
                                                                "#e7f5b6", 
                                                                "#3db2c4", 
                                                                "#0d3286",
                                                                "#043b5c"))(100)) +
              theme_feature() +
              facet_wrap(~Gene))
      
    } else {
      stop("Unknown reduction type")
    }
  }
}

Heatmap.Jo <- function(Seurat,
                      markers,
                      lev,
                      df,
                      num.genes = 5,
                      rowsize = 10,
                      annot) {
  
  if(missing(markers)) {
    
    markers = FindAllMarkers(Seurat, only.pos = T)
  }
  
  if(missing(lev)) {
    
    lev = levels(Seurat_DS@active.ident)
  }
  
  markers = markers %>%
    mutate(cluster = factor(cluster, levels = lev))
  
  df = df %>%
    mutate(Cluster = factor(Cluster, levels = lev))
  
  markers$pval.log = -log10(markers$p_val)
  markers$order = nrow(markers):1
  #obtain top20 differentially expressed genes for each cluster
  top10 = markers %>%
    group_by(cluster) %>%
    top_n(num.genes, order)
  top10$gene = as.character(top10$gene)
  
  #save heatmap
  require(Seurat)
  dat.heat = DoHeatmap(object = Seurat,
                       features = top10$gene)
  
  #create color annotations
  dat.heat.df = dat.heat$data[, -4]
  dat.heat.df = suppressWarnings(dcast(dat.heat.df, Feature~Cell, value.var = "Expression")) %>%
    column_to_rownames(var = "Feature")
  gene.ident = top10 %>%
    arrange(cluster, desc(pval.log)) %>%
    dplyr::select(c(gene, cluster, pval.log)) %>%
    filter(gene %in% rownames(dat.heat.df))
  gene.ident = gene.ident[!duplicated(gene.ident$gene),]
  cell.ident = df %>%
    dplyr::select(c(Cell, Cluster)) %>%
    arrange(Cluster)
  m = suppressWarnings(data.frame(gene = row.names(dat.heat.df)) %>%
    left_join(gene.ident, by = "gene") %>%
    arrange(cluster, desc(pval.log))) %>%
    .[complete.cases(.),]
  dat.heat.df = dat.heat.df[unique(m$gene), cell.ident$Cell]
  annotation_col = cell.ident %>%
    column_to_rownames(var = "Cell")
  
  if(missing(annot)) {
    

  
  #plot heatmap
  p = pheatmap::pheatmap(as.matrix(dat.heat.df),
               cluster_cols = F,
               cluster_rows = F,
               show_colnames = F,
               annotation_col = annotation_col,
               color = colorRampPalette(c("#cdcdcb", "#e7f5b6", "#3db2c4", "#0d3286"))(100),
               fontsize_row = rowsize,
               gaps_col = head(as.numeric(cumsum(table(df$Cluster))), -1),
               gaps_row = head(as.numeric(cumsum(table(gene.ident$cluster))), -1))
  
  p
  
  } else {
  
    ann_colors = list(Cluster = annot)
    p = pheatmap::pheatmap(as.matrix(dat.heat.df),
               cluster_cols = F,
               cluster_rows = F,
               show_colnames = F,
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               color = colorRampPalette(c("#cdcdcb", "#e7f5b6", "#3db2c4", "#0d3286"))(100),
               fontsize_row = rowsize,
               gaps_col = head(as.numeric(cumsum(table(df$Cluster))), -1),
               gaps_row = head(as.numeric(cumsum(table(gene.ident$cluster))), -1))
  }
  
  p
}

branch_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds) == FALSE]
  return(branch_points)
}

leaf_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds) == FALSE]
  return(leaves)
}

root_nodes <- function(cds, reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                           cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}

my.aggregate.Matrix = function (x, groupings = NULL, form = NULL, fun = "sum", ...)
{
  if (!methods::is(x, "Matrix"))
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  if (fun == "count")
    x <- x != 0
  groupings2 <- data.frame(A=as.factor(groupings))
  if (is.null(form))
    form <- stats::as.formula("~0+.")
  form <- stats::as.formula(form)
  mapping <- Matrix.utils::dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- Matrix::t(mapping) %*% x
  if (fun == "mean")
    result <- result/as.numeric(table(groupings)[rownames(result)])
  attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
                                                             groupings2$A))
  return(result)
}
