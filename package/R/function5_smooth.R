# functions for scClone
# bsh2024


#' @title Compute K value
#' @description Estimate signature-based K for a given variant.
#' @param var_string Variant string "chr_pos_ref>alt".
#' @param ref Genome build hg19/hg38/mm10.
#' @param sig Signature 96-vector.
#' @param K_smooth Smoothing parameter.
#' @return K value between 0 and 1.
#' @export
infer_K <- function(var_string, ref = ref, sig = sig, K_smooth = 0.2){
  var_array <- unlist(strsplit(var_string, split = "_"))
  var_array <- c(var_array[1:2],unlist(strsplit(var_array[3], ">"))[1:2])
  if(var_array[3] %in% c("A","T","C","G") & var_array[4] %in% c("A","T","C","G")){
    sca_vr = VariantAnnotation::VRanges(
      seqnames =  var_array[1],
      ranges = IRanges::IRanges(start = as.numeric(var_array[2]),end = as.numeric(var_array[2])),
      ref = var_array[3],
      alt = var_array[4])
    if(ref=="mm10"){
      sca_motifs = SomaticSignatures::mutationContext(sca_vr, BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
    }else if(ref=="hg19"){
      sca_motifs = SomaticSignatures::mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
    }else if(ref=="hg38"){
      sca_motifs = SomaticSignatures::mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    }else{
      cat("Error ref provided.\n")
    }
    var_96context <- as.data.frame(sca_motifs)
    var_96context$signature <- paste0(substring(var_96context$context,1,1),"[",substring(var_96context$alteration,1,1),">",
                                      substring(var_96context$alteration,2,2),"]",substring(var_96context$context,3,3))
    var_96context <- var_96context[,c(1:3,6,7,12:14)]
    K <- K_smooth+(1-K_smooth)*sig[var_96context$signature]/max(sig) #######################
  }else{
    K = 0.5
  }
  return(K)
}


#' @title Smooth VAF matrix
#' @description Smooth VAF values using local leaders/neighbors and signature weight.
#' @param matrix_VAF VAF matrix.
#' @param matrix_REF REF count matrix.
#' @param matrix_ALT ALT count matrix.
#' @param meta Metadata with celltype and cluster_tem.
#' @param leader_min_depth Min depth for leader cells.
#' @param neighbor_min_dis Max distance for neighbors.
#' @param ref Genome build.
#' @param sig Optional 96-signature vector.
#' @param method Distance metric.
#' @param K_smooth Smoothing parameter.
#' @return Smoothed VAF matrix.
#' @export
make_smooth_matrix <- function(matrix_VAF, matrix_REF, matrix_ALT, meta, leader_min_depth, neighbor_min_dis, ref, sig = NULL, method = method, K_smooth = 0.2){
  matrix_smooth <- matrix_VAF
  dis <- make_distance_matrix(matrix_VAF = matrix_VAF, method = method)
  pb <- txtProgressBar(min = 0, max = nrow(matrix_VAF), style = 3)
  options(warn = -1)
  flag = 0
  for(i in 1:nrow(matrix_VAF)){
    setTxtProgressBar(pb, i)
    cat(paste0(" \r",flag," VAFs smoothed"),"")
    var_string <- rownames(matrix_VAF)[i]
    if(!is.null(sig)==T){
      K <- infer_K(var_string, ref = ref, sig = sig, K_smooth = K_smooth)
    }else{
      K = 0.5
    }

    # set leader
    depth <- as.numeric(matrix_ALT[i,]) + as.numeric(matrix_REF[i,])
    leaders <- colnames(matrix_VAF)[which(depth>=leader_min_depth)]
    leaders_cell_type <- as.character(meta$celltype[match(leaders, meta$cell)])
    leaders_cluster <- as.character(meta$cluster_tem[match(leaders, meta$cell)])
    for(j in 1:ncol(matrix_VAF)){
      # set cell
      cell <- colnames(matrix_VAF)[j]
      cell_type <- as.character(meta$celltype[which(meta$cell==cell)])
      cluster <- as.character(meta$cluster_tem[which(meta$cell==cell)])
      # set neighbor
      neighbors <- rownames(dis)[which( dis[,cell]>0 & dis[,cell]<neighbor_min_dis )]
      neighbors_cell_type <- meta$celltype[match(neighbors, meta$cell)]
      neighbors_cluster <- meta$cluster_tem[match(neighbors, meta$cell)]
      # final leader & neighbor
      if(length(cluster)==0){
        neighbors <- neighbors[which(neighbors_cell_type==cell_type)]
        leaders_tem <- leaders[which(leaders_cell_type==cell_type)]
      }else{
        neighbors <- neighbors[which(neighbors_cell_type==cell_type & neighbors_cluster==cluster)]
        leaders_tem <- leaders[which(leaders_cell_type==cell_type & leaders_cluster==cluster)]
      }
      # smooth
      if(length(intersect(neighbors,leaders_tem))>0){
        # ÓÐleader
        leaders_tem <- intersect(neighbors,leaders_tem)
        smooth <- mean(as.numeric(matrix_VAF[i,which(colnames(matrix_VAF) %in% leaders_tem)]))
        matrix_smooth[i,j] <- smooth
        flag = flag + 1
      }else if(length(neighbors)>0){
        # ÎÞleader ÓÃneighbor
        smooth <- as.numeric(matrix_VAF[i,which(colnames(matrix_VAF) %in% neighbors)])
        smooth <- matrix_VAF[i,j]*(1-K) + mean(smooth)*K
        matrix_smooth[i,j] <- smooth
        flag = flag + 1
      }
    }
  }
  close(pb)
  cat("DONE! \n")
  return(matrix_smooth)
}
