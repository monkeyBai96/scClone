# 10000 vcf.gz
# bsh 0515
rm(list = ls())
graphics.off()
source("H:/scclone/function0_color.R")
source("H:/scclone/function0_statistic.R")
source("H:/scclone/function1_read_data.R")
source("H:/scclone/function2_signature.R")
source("H:/scclone/function3_make_vaf.R")
source("H:/scclone/function4_cluster.R")
source("H:/scclone/function5_smooth.R")
source("H:/scclone/function6_randomforest.R")
sample2 <- c("P4_1","P6_1")
cutoff_perc <- c(0.9,0.9)
# color
celltypes <- paste0("type_",c(0:11))
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(8))
show_col(celltypes_color)
celltypes_color <- celltypes_color[c(1:5,13,12,9,11,6,7,10)]
names(celltypes_color) <- celltypes
for(i in 1:2){
  ################################## path ########################################
  path = "H:/scclone/PMID32579974_ST"
  rna_path = paste0(path,"/paper/ST/GSE144239_RAW/",sample2[i])
  mut_path = paste0(path,"/",sample2[i],"/annovar")
  work_path = paste0(path,"/",sample2[i],"/work")
  setwd(paste0(work_path,"/5_smooth"))
  ################################################################################
  
  # read
  meta <- read.csv(paste0(work_path, "/1_basic/meta_select.csv"), row.names = 1)
  matrix_VAF <- read.csv(paste0(work_path, "/4_vaf/matrix_VAF.csv"), row.names = 1)
  matrix_ALT <- read.csv(paste0(work_path, "/4_vaf/matrix_ALT.csv"), row.names = 1)
  matrix_REF <- read.csv(paste0(work_path, "/4_vaf/matrix_REF.csv"), row.names = 1)
  matrix_VAF[is.na(matrix_VAF)] <- 0
  matrix_VAF[matrix_VAF>2] <- 2
  hist(unlist(matrix_VAF), xlim = c(0.1,3), breaks = 1000)
  
  # hclust
  demo <- as.matrix(dist(t(matrix_VAF), method = "manhattan", diag = FALSE, upper = FALSE))
  hist(unlist(demo), xlim = c(0,max(demo)), breaks = 100)
  #cutoff <- quantile(unlist(demo),cutoff_perc[i])
  cluster <- get_hclust(matrix = matrix_VAF, n_cluster = 2, method = "manhattan")
  meta$cluster_tem <- paste0("cluster_",cluster[match(meta$cell, names(cluster))])
  table(meta$cluster_tem)
  
  # refine meta: celltype and cluster_tem
  meta_order <- meta %>% arrange(celltype)
  meta_select <- refine_meta(meta = meta_order, min_cell = 5)
  table(meta_select$cluster_tem, meta_select$celltype_resolution)
  
  # color
  anno_color <- list(celltype = celltypes_color)
  p1 <- pheatmap(matrix_VAF, annotation_col = meta_order[,c("cluster_tem","celltype","celltype_resolution")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
  p2 <- pheatmap(matrix_VAF[,meta_order$cell], annotation_col = meta_order[,c("cluster_tem","celltype","celltype_resolution")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = F)
  
  # signature
  signature <- read.csv(paste0(work_path,"/2_svm/signature.csv"))
  sig <- table(signature$signature)
  
  # slow
  matrix_VAF <- matrix_VAF[,which(colnames(matrix_VAF) %in% rownames(meta_select))]
  matrix_ALT <- matrix_ALT[,which(colnames(matrix_ALT) %in% rownames(meta_select))]
  matrix_REF <- matrix_REF[,which(colnames(matrix_REF) %in% rownames(meta_select))]
  meta_select$celltype <- meta_select$celltype_resolution ##################################
  matrix_smooth <- make_smooth_matrix(matrix_VAF = matrix_VAF, matrix_ALT = matrix_ALT,matrix_REF = matrix_REF,meta = meta_select, sig = sig,
                                      ref = "hg38", 
                                      method = "manhattan",
                                      leader_min_depth = 10,
                                      neighbor_min_dis = quantile(unlist(demo),0.5),
                                      K_smooth = 0.4)
  graphics.off()
  meta_select$celltype <- paste0("type_", meta_select$Spatial_snn_res.0.8)
  p3 <- pheatmap(matrix_smooth, annotation_col = meta_select[,c("cluster_tem","celltype","celltype_resolution")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
  p4 <- pheatmap(matrix_smooth[,meta_select$cell], annotation_col = meta_select[,c("cluster_tem","celltype","celltype_resolution")], colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = F)
  
  # save
  setwd(paste0(work_path,"/5_smooth"))
  write.csv(matrix_smooth, "matrix_smooth.csv", row.names = T)
  write.csv(matrix_VAF, "matrix_VAF.csv", row.names = T)
  write.csv(matrix_ALT, "matrix_ALT.csv", row.names = T)
  write.csv(matrix_REF, "matrix_REF.csv", row.names = T)
  write.csv(meta_select, paste0(work_path,"/1_basic/meta_select.csv"), row.names = T)
  p <- plot_grid(p1$gtable, p2$gtable, p3$gtable, p4$gtable, ncol = 2)
  ggsave("heatmap_smooth.pdf", p, width = 20, height = 20)
}














