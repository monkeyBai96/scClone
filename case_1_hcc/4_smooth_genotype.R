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
HCC_name <- c("HCC1","HCC2","HCC5","HCC9")
cutoff_perc <- c(0.9,0.8,0.9,0.9)
for(i in 2:3){
  ################################## path ########################################
  path = "H:/scclone/PMID_SU_HCC"
  rna_path = paste0(path,"/paper")
  mut_path = paste0(path,"/",HCC_name[i],"/annovar")
  work_path = paste0(path,"/",HCC_name[i],"/work")
  setwd(paste0(work_path,"/5_smooth"))
  ################################################################################
  
  
  # read data
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  matrix_VAF <- read.csv(paste0(work_path,"/4_vaf/matrix_VAF.csv"), row.names = 1)
  matrix_ALT <- read.csv(paste0(work_path,"/4_vaf/matrix_ALT.csv"), row.names = 1)
  matrix_REF <- read.csv(paste0(work_path,"/4_vaf/matrix_REF.csv"), row.names = 1)
  matrix_VAF[is.na(matrix_VAF)] <- 0
  matrix_VAF[matrix_VAF>2] <- 2
  
  # hclust
  demo <- as.matrix(dist(t(matrix_VAF), method = "manhattan", diag = FALSE, upper = FALSE))
  hist(unlist(demo), xlim = c(0,max(demo)), breaks = 100)
  cutoff <- quantile(unlist(demo),cutoff_perc[i])
  cluster <- get_hclust(matrix = matrix_VAF, cut = cutoff, method = "manhattan")
  meta$cluster_tem <- paste0("cluster_",cluster[match(meta$cell, names(cluster))])
  table(meta$cluster_tem)
  
  # refine meta: celltype and cluster_tem
  meta_order <- meta %>% arrange(celltype)
  meta_select <- refine_meta(meta = meta_order, min_cell = 5)
  table(meta_select$cluster_tem)
  
  # color
  n = length(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell"))
  color <- c(pal_npg("nrc")(10))[c(1:n)]
  show_col(color)
  names(color) <- c(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell"))
  anno_color <- list(celltype = color)
  p1 <- pheatmap(matrix_VAF, annotation_col = meta_order[,c("cluster_tem","celltype","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
  p2 <- pheatmap(matrix_VAF[,meta_order$cell], annotation_col = meta_order[,c("cluster_tem","celltype","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = F)
  
  # signature
  signature <- read.csv(paste0(work_path,"/2_svm/signature.csv"))
  sig <- table(signature$signature)
  
  # slow
  matrix_VAF <- matrix_VAF[,which(colnames(matrix_VAF) %in% rownames(meta_select))]
  matrix_ALT <- matrix_ALT[,which(colnames(matrix_ALT) %in% rownames(meta_select))]
  matrix_REF <- matrix_REF[,which(colnames(matrix_REF) %in% rownames(meta_select))]
  #meta_select$cluster_tem <- meta_select$sample ##################################
  matrix_smooth <- make_smooth_matrix(matrix_VAF = matrix_VAF, matrix_ALT = matrix_ALT,matrix_REF = matrix_REF,meta = meta_select, sig = sig,
                                      ref = "hg19", 
                                      method = "manhattan",
                                      leader_min_depth = 50,
                                      neighbor_min_dis = quantile(unlist(demo),0.5),
                                      K_smooth = 0.2)
  graphics.off()
  p3 <- pheatmap(matrix_smooth, annotation_col = meta_select[,c("cluster_tem","celltype","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
  p4 <- pheatmap(matrix_smooth[,meta_select$cell], annotation_col = meta_select[,c("cluster_tem","celltype","sample")], colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
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







# plot VAF calculate
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
HCC_name <- c("HCC1","HCC2","HCC5","HCC9")
cutoff_perc <- c(0.9,0.8,0.9,0.9)
for(i in 2:3){
  ################################## path ########################################
  path = "H:/scclone/PMID_SU_HCC"
  rna_path = paste0(path,"/paper")
  mut_path = paste0(path,"/",HCC_name[i],"/annovar")
  work_path = paste0(path,"/",HCC_name[i],"/work")
  setwd(paste0(work_path,"/5_smooth"))
  ################################################################################
  
  matrix_ALT <- read.csv("matrix_ALT.csv", row.names = 1)
  matrix_REF <- read.csv("matrix_REF.csv", row.names = 1)
  matrix_smooth <- read.csv("matrix_smooth.csv", row.names = 1)
  single <- read.csv(paste0(work_path, "/3_filter/single_filter_more.csv"))
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  
  single$VAF <- single$AD_2/(single$AD_1+single$AD_2)
  matrix_VAF_raw <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","VAF")], var~cell, value.var = "VAF", fun.aggregate=sum)))
  matrix_VAF_raw[is.na(matrix_VAF_raw)] <- 0
  matrix_VAF_raw <- matrix_VAF_raw[rownames(matrix_smooth),colnames(matrix_smooth)]
  data_tem <- data.frame(VAF = as.numeric(unlist(matrix_VAF_raw)),
                         ALT = as.numeric(unlist(matrix_ALT)),
                         REF = as.numeric(unlist(matrix_REF)),
                         genotype = as.numeric(unlist(matrix_smooth)))
  data_tem <- data_tem %>% mutate(alt_depth = ifelse(ALT==0, "drop-out",ifelse(ALT+REF<10,"<10",">10")))
  data_tem$alt_depth <- factor(data_tem$alt_depth, levels = c(">10","<10","drop-out"))
  color <- pal_npg("nrc")(10)[c(2,1,3)]
  p1 <- ggplot(data_tem, aes(VAF, fill = alt_depth, group = alt_depth))+
    geom_histogram(position = "stack")+
    scale_fill_manual(values = color)+
    theme_few()
  p2 <- ggplot(data_tem, aes(genotype, fill = alt_depth, group = alt_depth))+
    geom_histogram(position = "stack")+
    scale_fill_manual(values = color)+
    theme_few()
  graphics.off()
  p <- gg.gap(plot=p1,segments=list(c(10000,20000)),ylim=c(0,3*10^5))
  #add.legend(plot = p1,margin = c(top=1,right=1,bottom=200,left=300))
  ggsave("hist_VAF_raw.pdf",p, width = 6, height = 5)
  
  p <- gg.gap(plot=p2,segments=list(c(30000,40000)),ylim=c(0,3*10^5))
  #add.legend(plot = p2,margin = c(top=1,right=1,bottom=200,left=300))
  ggsave("hist_VAF.pdf",p, width = 6, height = 5)
  
  # plot distance
  meta_order <- meta %>% arrange(celltype)
  matrix_order <- matrix_smooth[,rownames(meta_order)]
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
  # color
  n = length(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell"))
  color <- c(pal_npg("nrc")(10))[c(1:n)]
  show_col(color)
  names(color) <- c(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell"))
  anno_color <- list(celltype = color)
  p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("cluster_tem","celltype","sample")], annotation_colors = anno_color,
                          annotation_row = meta_order[,c("cluster_tem","celltype","sample")],
                          show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                          color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
  ggsave("distance.pdf", p, width = 10, height = 10)

}







