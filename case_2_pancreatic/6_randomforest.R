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
sample2 <- c("P3_YF","P3_ZY")
cutoff_perc <- c(0.9,0.9)
# color
celltypes <- c("B_cell","Ductal_cell","Fibroblasts","Mac/Mono","Mast_cell","MKI67+cell","NK","Plasma_cell","T_cell")
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(1:9)]
names(celltypes_color) <- celltypes
for(i in 1:1){
  ################################## path ########################################
  path = "H:/scclone/PMID37612267"
  rna_path = paste0(path,"/paper")
  mut_path = paste0(path,"/",sample2[i],"/annovar")
  work_path = paste0(path,"/",sample2[i],"/work")
  setwd(paste0(work_path,"/7_rf"))
  #################################################################################
  
  
  # read data
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  matrix_VAF <- read.csv(paste0(work_path, "/6_rpca/matrix_rpca.csv"), row.names = 1)
  matrix_VAF <- rounding_VAF(matrix_VAF)
  
  # run randomforest
  demo <- run_randomforest(matrix = matrix_VAF, meta = meta, train_ratio = 0.7, output_path = paste0(work_path,"/7_rf"))
  saveRDS(demo,"RF_result.rds")
  importance <- demo[["importance"]]
  var_keep <- rownames(importance)[which(importance$MeanDecreaseAccuracy>0&importance$MeanDecreaseGini>0)]
  matrix_rf <- matrix_VAF[var_keep,]
  meta_order <- meta %>% arrange(celltype)
  matrix_rf_order <- matrix_rf[,meta_order$cell]
  
  # color
  anno_color <- list(celltype = celltypes_color)
  graphics.off()
  p1 <- pheatmap(matrix_rf, annotation_col = meta[,c("cluster_tem","celltype","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50),annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
  p2 <- pheatmap(matrix_rf_order, annotation_col = meta_order[,c("cluster_tem","celltype","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50),annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = F)
  p <- plot_grid(p1$gtable, p2$gtable, ncol = 2)
  ggsave("heatmap_rf.pdf",p, width = 20, height = 10)
  write.csv(matrix_rf, "matrix_rf.csv", row.names = T)
  
  # plot distance
  meta_order <- meta %>% arrange(celltype)
  matrix_order <- matrix_rf[,rownames(meta_order)]
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
  p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("cluster_tem","celltype","sample")], annotation_colors = anno_color,
                          annotation_row = meta_order[,c("cluster_tem","celltype","sample")],
                          show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                          color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
  ggsave("distance.pdf",p, width = 10, height = 10)
}
