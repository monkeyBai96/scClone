# functions for scClone
# bsh2024


#' @title Compute distance matrix
#' @description Calculate cell-to-cell distance from VAF matrix.
#' @param matrix_VAF VAF matrix (variants × cells).
#' @param method Distance metric.
#' @return Symmetric distance matrix.
#' @export
make_distance_matrix <- function(matrix_VAF = matrix_VAF, method = method){
  dis <- as.matrix(dist(t(matrix_VAF), method = method, diag = FALSE, upper = FALSE))
  return(dis)
}


#' @title Hierarchical clustering
#' @description Cluster cells by VAF distance and cut tree.
#' @param matrix VAF matrix.
#' @param method Distance metric.
#' @param n_cluster Number of clusters.
#' @param cut Height cutoff.
#' @return Vector of cluster labels.
#' @export
get_hclust <- function(matrix = matrix, method = method, n_cluster = NULL, cut = NULL){
  dist <- suppressMessages(dist(t(matrix), method = method, diag = FALSE, upper = FALSE))
  model <- hclust(dist, method = "complete", members = NULL)
  cluster <- cutree(model, k = n_cluster, h = cut)
  return(cluster)
}



#' @title Celltype-wise clustering
#' @description Cluster cells within each celltype and append new cluster labels.
#' @param matrix VAF matrix.
#' @param meta Metadata with celltype.
#' @return Metadata with appended cluster_new column.
#' @export
get_hclust_celltype <- function(matrix = matrix, meta = meta){
  matrix <- matrix[rowSums(matrix)>0,]
  new_cluster <- c()
  for(i in unique(meta$celltype)){
    cell_select <- meta$cell[which(meta$celltype==i)]
    tem <- t(matrix[,cell_select])
    capture.output(result <- try({
      eclust(tem, "hclust",k.max=ceiling(log(length(cell_select), base=5)), hc_metric = "manhattan", hc_method = "complete", stand = T)
    }, silent = TRUE))

    if(inherits(result, "try-error")){
      tem_cluster <- rep(1,length(cell_select))
      names(tem_cluster) <- rownames(tem)
      new_cluster <- c(new_cluster,tem_cluster)
    }else if(length(unique(result$cluster))==1){
      tem_cluster <- rep(1,length(cell_select))
      names(tem_cluster) <- rownames(tem)
      new_cluster <- c(new_cluster,tem_cluster)
    }else{
      tem_cluster <- result$cluster
      count <- data.frame(cell = names(tem_cluster), cluster = as.character(tem_cluster), count = as.numeric(colSums(matrix)[names(tem_cluster)]))
      count_cluster <- aggregate(count ~ cluster, data = count, FUN = mean) %>% arrange(count)
      count$cluster_order <- match(tem_cluster,count_cluster$cluster)
      tem_cluster <- count$cluster_order
      names(tem_cluster) <- count$cell
      new_cluster <- c(new_cluster,tem_cluster)
    }
    cat(paste0(i,": ", length(cell_select)," cells; ",max(tem_cluster)," clusters.\n"))
  }
  meta <- meta[names(new_cluster),]
  meta$cluster_new <- paste0(meta$celltype,"_",new_cluster)
  return(meta)
}


#' @title Refine metadata
#' @description Remove low-count celltypes/clusters and re-label ambiguous clusters.
#' @param meta Metadata data frame.
#' @param min_cell Minimum cell count per group.
#' @return Updated metadata with cluster_tem refined.
#' @export
refine_meta <- function(meta = meta, min_cell = min_cell){
  table_celltype <- table(meta$celltype)
  keep_celltype <- names(table_celltype)[which(table_celltype > min_cell)]
  meta <- meta %>% filter(celltype %in% keep_celltype)
  table_cluster <- table(meta$cluster_tem)
  keep_cluster <- names(table_cluster)[which(table_cluster > min_cell)]
  meta$cluster_tem <- ifelse(meta$cluster_tem %in% keep_cluster, meta$cluster_tem, NA)

  table <- as.data.frame(table(meta$celltype, meta$cluster_tem))
  table <- dcast(table, value.var = "Freq", formula = Var1~Var2)
  table$Sum <- table(meta$celltype)[match(table$Var1, names(table(meta$celltype)))]
  table <- col1_to_rowname(table)
  table <- rbind(table,
                 simpson = vegan::diversity(t(table), index = "simpson"),
                 shannon = vegan::diversity(t(table), index = "shannon"))
  keep_cluster <- unique(colnames(table)[which(table["simpson",-ncol(table)] < 1.2*table["simpson",ncol(table)])],
                         colnames(table)[which(table["shannon",-ncol(table)] < 1.2*table["shannon",ncol(table)])])
  meta$cluster_tem <- ifelse(meta$cluster_tem %in% keep_cluster, meta$cluster_tem, NA)
  return(meta)
}


#' @title Finalize cluster labels
#' @description Re-label clusters within each celltype sequentially.
#' @param meta Metadata data frame.
#' @return Metadata with sequential cluster_new labels.
#' @export
refine_meta_final <- function(meta = meta){
  list <- list()
  celltype = unique(meta$celltype)
  for(c in 1:length(celltype)){
    tem <- meta[which(meta$celltype==celltype[c]),]
    if(length(unique(tem$cluster_new))>1){
      n <- length(unique(tem$cluster_new))
      tem <- tem %>% arrange(cluster_new)
      cluster_new <- unique(tem$cluster_new)
      for(cl in 1:length(cluster_new)){
        tem$cluster_new[which(tem$cluster_new==cluster_new[cl])] <- paste0(celltype[c],"_",cl)
      }
    }else{
      tem$cluster_new <- paste0(tem$celltype,"_1")
    }
    list[[c]] <- tem
  }
  meta <- do.call(rbind,list)
  return(meta)
}



