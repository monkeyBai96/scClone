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
  setwd(paste0(work_path,"/6_rpca"))
  ################################################################################
  
  
  ## read data
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  matrix_VAF <- read.csv(paste0(work_path,"/5_smooth/matrix_smooth.csv"), row.names = 1)
  res <- rpca(as.matrix(matrix_VAF), max.iter = 500)
  saveRDS(res, "rpca.rds")
  
  res <- readRDS("rpca.rds")
  colnames(res$L) <- colnames(matrix_VAF)
  colnames(res$S) <- colnames(matrix_VAF)
  # color
  anno_color <- list(celltype = celltypes_color)
  p1 <- pheatmap(matrix_VAF, show_rownames = F, show_colnames = F, annotation_col = meta[,c("cluster_tem","celltype","celltype_resolution")],color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color)
  p2 <- pheatmap(res$L, show_rownames = F, show_colnames = F, annotation_col = meta[,c("cluster_tem","celltype","celltype_resolution")],color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color)
  p3 <- pheatmap(res$S, show_rownames = F, show_colnames = F, annotation_col = meta[,c("cluster_tem","celltype","celltype_resolution")],color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color)
  p <- plot_grid(p1$gtable, p2$gtable, p3$gtable, ncol = 3)
  ggsave("heatmap_rpca.pdf",p, width = 30, height = 10)
  matrix_rpca <- res$L
  rownames(matrix_rpca) <- rownames(matrix_VAF)
  write.csv(matrix_rpca, "matrix_rpca.csv")
  
  # plot distance
  meta_order <- meta %>% arrange(celltype)
  matrix_order <- matrix_rpca[,rownames(meta_order)]
  distance <- as.matrix(suppressMessages(dist(t(matrix_order), method = "manhattan", diag = FALSE, upper = FALSE)))
  order <- c()
  group <- unique(meta_order$celltype)
  for(j in 1:length(group)){
    d <- dist(t(matrix_order[,which(meta_order$celltype==group[j])]), method = "manhattan", diag = FALSE, upper = FALSE)
    if(nrow(d)>2){
      x <- hclust(d, method = "complete")
      order <- c(order,x$labels[x$order])
    }else{
      order <- c(order,meta_order$cell[which(meta_order$celltype==group[j])])
    }
  }
  p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("cluster_tem","celltype","celltype_resolution")], annotation_colors = anno_color,
                          annotation_row = meta_order[,c("cluster_tem","celltype","celltype_resolution")],
                          show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                          color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
  ggsave("distance.pdf",p, width = 10, height = 10)
}

