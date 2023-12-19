celltype.heatmap.local <- function (data, assay = "RNA", genes, ref = NULL, scale = "row", 
                                    method = c("ward.D2", "ward.D", "average"), brewer.palette = "RdBu", 
                                    palette_reverse = F, palette = NULL, cluster.col = "functional.cluster", 
                                    metadata = NULL, order.by = NULL, flip = FALSE, cluster_genes = FALSE, 
                                    cluster_samples = FALSE, min.cells = 10, show_samplenames = FALSE,
                                    return.matrix = FALSE,
                                    remove.NA.meta = TRUE, breaks = seq(-2, 2, by = 0.1), ...) 
{
  set.seed(123)
  method = method[1]
  if (is.null(metadata)) {
    stop("Must at least provide one metadata")
  }
  else {
    meta.sub <- data@meta.data[, which(colnames(data@meta.data) %in% 
                                         cluster.col), drop = F]
    meta.sub <- cbind(meta.sub, data@meta.data[, metadata, 
                                               drop = F])
  }
  meta.sub[meta.sub == "NA"] = NA
  if (remove.NA.meta == TRUE) {
    meta.sub <- drop_na(meta.sub)
  }
  data$metaSubset <- factor(apply(meta.sub, 1, paste, collapse = "!"))
  t <- table(data$metaSubset)
  accept <- names(t)[t > min.cells]
  data <- subset(data, subset = metaSubset %in% accept)
  m <- c()
  genes.removed <- c()
  for (g in unique(genes)) {
    if (g %in% rownames(data@assays[[assay]])) {
      m[[g]] <- tapply(data@assays[[assay]][g, ], data$metaSubset, 
                       mean)
    }
    else {
      genes.removed <- c(genes.removed, g)
    }
  }
  if (length(genes.removed) > 0) {
    cat("These genes were not found in the assay and were excluded from plotting:", 
        genes.removed)
  }
  m <- as.data.frame(m)
  m <- m[accept, ]
  m.subset <- factor(unlist(lapply(strsplit(rownames(m), "!", 
                                            perl = T), function(x) x[[1]])))
  m.meta <- list()
  for (i in 1:length(metadata)) {
    m.meta[[i]] <- factor(unlist(lapply(strsplit(rownames(m), 
                                                 "!", perl = T), function(x) x[[i + 1]])))
  }
  names(m.meta) <- metadata
  m.meta <- as.data.frame(m.meta)
  m <- cbind(m, m.subset, m.meta)
  if (!is.null(ref)) {
    m$m.subset <- factor(m$m.subset, levels = levels(ref$functional.cluster))
    m <- arrange(m, m.subset)
    m.subset <- factor(unlist(lapply(strsplit(rownames(m), 
                                              "!", perl = T), function(x) x[[1]])))
  }
  if (!is.null(order.by)) {
    m <- arrange(m, m.subset, get(order.by))
    m.meta <- list()
    for (i in 1:length(metadata)) {
      m.meta[[i]] <- factor(unlist(lapply(strsplit(rownames(m), 
                                                   "!", perl = T), function(x) x[[i + 1]])))
    }
    names(m.meta) <- metadata
  }
  m <- m[1:(length(m) - length(metadata) - 1)]
  require(RColorBrewer)
  if (palette_reverse) {
    color = colorRampPalette(rev(brewer.pal(n = 7, name = brewer.palette)))(length(breaks))
  }
  else {
    color = colorRampPalette(brewer.pal(n = 7, name = brewer.palette))(length(breaks))
  }
  palettes.default <- c("Paired", "Set2", "Accent", "Dark2", 
                        "Set1", "Set3")
  if (is.null(palette)) {
    palette <- list()
    palette[["Subtype"]] <- colorRampPalette(brewer.pal(n = 8, 
                                                        name = "Set1"))(length(unique(unlist(m.subset))))
    names(palette[["Subtype"]]) <- c(unique(m.subset))
    for (i in 1:length(metadata)) {
      meta <- metadata[i]
      palette[[meta]] <- colorRampPalette(brewer.pal(n = 6, 
                                                     name = palettes.default[i]))(length(unique(m.meta[[meta]])))
      names(palette[[meta]]) <- levels(m.meta[[meta]])
    }
  }
  if (!is.null(ref)) {
    palette[["Subtype"]] <- ref@misc$atlas.palette
  }
  annotation_col = data.frame(Subtype = m.subset)
  annotation_col <- cbind(annotation_col, m.meta)
  rownames(annotation_col) = rownames(m)
  
  if (return.matrix) {
    m <- t(m)
    mean <- apply(m, 1, mean, na.rm = T)
    sd <- apply(m, 1, sd, na.rm = T)
    return((m - mean) / sd)
  }
  
  if (flip) {
    h <- pheatmap::pheatmap(m, cluster_rows = cluster_samples, 
                            cluster_cols = cluster_genes, scale = scale, breaks = breaks, 
                            color = color, annotation_row = annotation_col, show_rownames = show_samplenames, 
                            border_color = NA, annotation_colors = palette, fontsize_row = 6, 
                            fontsize = 7, clustering_method = method, ...)
  }
  else {
    h <- pheatmap::pheatmap(t(m), cluster_rows = cluster_genes, 
                            cluster_cols = cluster_samples, scale = scale, breaks = breaks, 
                            color = color, annotation_col = annotation_col, show_colnames = show_samplenames, 
                            border_color = NA, annotation_colors = palette, fontsize_row = 6, 
                            fontsize = 7, clustering_method = method, ...)
  }
  return(h)
}


#Make a subset vs. study heatmap of pseudobulk profiles

make.heatmap <- function(data, assay="RNA", genes, scale="row",
                         method=c("ward.D2","ward.D", "average"), brewer.palette="RdBu",
                         cluster.col = "functional.cluster", dataset.col = "study", dataset.tissue= "tissue", flip=FALSE,
                         cluster_genes = FALSE, cluster_samples=FALSE, min.cells = 10,
                         palette = NULL) {
  
  require(pheatmap)
  set.seed(123)
  
  method = method[1]
  
  data$Study_Subset <- factor(paste(data@meta.data[,cluster.col],data@meta.data[,dataset.col],data@meta.data[,dataset.tissue],sep=","))
  
  t <- table(data$Study_Subset)
  accept <- names(t)[t>min.cells]
  
  data <- subset(data, subset=Study_Subset %in% accept)
  
  #Calculate mean expression by cluster
  m <- c()
  for( g in unique(genes)){
    m[[g]] <- tapply(data@assays[[assay]][g,],data$Study_Subset, mean)
  }
  m <- as.data.frame(m)
  
  m <- m[accept,]
  
  m.study <- factor(unlist(lapply(strsplit(rownames(m),"\\,",perl = T),function(x) x[[2]])))
  m.subset <- factor(unlist(lapply(strsplit(rownames(m),"\\,",perl = T),function(x) x[[1]])))
  m.tissue <- factor(unlist(lapply(strsplit(rownames(m),"\\,",perl = T),function(x) x[[3]])))
  #gaps_col <- (which(!duplicated(m.subset))-1)[-1]
  
  breaksList = seq(-2, 2, by = 0.1)
  
  require(RColorBrewer)
  color = colorRampPalette(rev(brewer.pal(n = 7, name = brewer.palette)))(length(breaksList))
  
  if (is.null(palette)) {
    palette = brewer.pal(n=length(unique(m.subset)), name="Paired")
    names(palette) <- unique(m.subset)
  }
  
  annotation_col = data.frame(
    Subtype = m.subset,
    Study = m.study,
    Tissue = m.tissue)
  rownames(annotation_col) = rownames(m)
  if (flip) { 
    h <- pheatmap::pheatmap(m, cluster_rows = cluster_samples,
                            cluster_cols = cluster_genes,scale = scale,
                            breaks = breaksList, color=color, 
                            annotation_row = annotation_col, 
                            show_rownames = F,
                            show_colnames = T,
                            border_color = NA,
                            annotation_colors = list(Subtype=palette), 
                            fontsize_row=6,fontsize = 7, 
                            clustering_method=method)
  } else {
    h <- pheatmap::pheatmap(t(m),cluster_rows = cluster_genes,
                            cluster_cols = cluster_samples,scale = scale,
                            breaks = breaksList, color=color, 
                            annotation_col = annotation_col, 
                            show_colnames = F,
                            show_rownames = T,
                            border_color = NA,
                            annotation_colors = list(Subtype=palette), 
                            fontsize_row=6,fontsize = 7, 
                            clustering_method=method)
  }
  return(h)
}

