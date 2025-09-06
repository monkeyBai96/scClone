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
sample7 <- c("P1","P2","P3","P4","P5","P6","P8")
for(i in c(3,7)){
  ################################################################################
  path = "H:/scclone/PMID38570491_ST"
  rna_path = paste0(path,"/paper/ST/",sample7[i])
  mut_path = paste0(path,"/",sample7[i],"/annovar")
  work_path = paste0(path,"/",sample7[i],"/work")
  setwd(paste0(work_path,"/3_filter"))
  ################################################################################
  
  
  # read data
  meta <- read.csv(paste0(work_path, "/0_exp/meta.csv"), row.names = 1)
  single_raw <- read.csv(paste0(work_path, "/1_basic/single_mut_vaf.csv"))
  single_raw <- single_raw %>% filter(cell %in% meta$file)
  single_raw$sample <- meta$sample[match(single_raw$cell, meta$file)]
  single_raw$cell <- meta$cell[match(single_raw$cell, meta$file)]
  single_raw$var_pos <- paste0(single_raw$Chr,"_",single_raw$Pos)
  single_raw$var <- paste0(single_raw$Chr,"_",single_raw$Pos,"_",single_raw$Ref,">",single_raw$Alt)
  
  single_svm <- read.csv(paste0(work_path, "/2_svm/single_svm.csv"))
  single_svm <- single_svm %>% filter(mark=="True" | mark1=="True")
  single_svm$var_pos <- paste0(single_svm$Chr,"_",single_svm$Pos)
  
  single_raw <- single_raw %>% filter(var_pos %in% single_svm$var_pos) %>% filter(Chr %in% paste0("chr",c(1:22)))
  single_raw$DP <- ifelse(!is.na(single_raw$DP),single_raw$DP,single_raw$DPI)
  single_raw <- split_VAF_col(single_raw, col = "AD")
  single_raw <- single_raw[,c(1:9,match(c("DP","AD_1","AD_2","cell","sample","var_pos","var"), colnames(single_raw)))]
  
  single_deep <- single_raw %>% filter(DP>10)
  
  # filter
  single <- mutation_filter(mut = single_deep, min_celltype = 8, min_mut_per_cell = 8, min_cell_one_mut = 5, meta = meta)
  meta_select <- meta %>% filter(cell %in% single$cell)
  table(meta_select$sample, meta_select$celltype)
  
  # save
  write.csv(single, "single_filter.csv", row.names = F)
  single_raw <- single_raw %>% 
    filter(var %in% single$var) %>%
    filter(cell %in% single$cell)
  write.csv(single_raw, "single_filter_more.csv", row.names = F)
  write.csv(meta_select, paste0(work_path,"/1_basic/meta_select.csv"), row.names = T)
}






# after filter
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
sample7 <- c("P1","P2","P3","P4","P5","P6","P8")
# color
celltypes <- paste0("type_",c(0:11))
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(8))
show_col(celltypes_color)
celltypes_color <- celltypes_color[c(1:5,13,12,9,11,6,7,10)]
names(celltypes_color) <- celltypes
for(i in c(3,7)){
  ################################################################################
  path = "H:/scclone/PMID38570491_ST"
  rna_path = paste0(path,"/paper/ST/",sample7[i])
  mut_path = paste0(path,"/",sample7[i],"/annovar")
  work_path = paste0(path,"/",sample7[i],"/work")
  setwd(paste0(work_path,"/4_vaf"))
  ################################################################################
  
  
  single <- read.csv(paste0(work_path, "/3_filter/single_filter_more.csv"))
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  
  # make matrix
  single$AD_1 <- as.numeric(single$AD_1)
  matrix_REF <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","AD_1")], var~cell, value.var = "AD_1", fun.aggregate=sum)))
  single$AD_2 <- as.numeric(single$AD_2)
  matrix_ALT <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","AD_2")], var~cell, value.var = "AD_2", fun.aggregate=sum)))
  dim(matrix_REF)
  dim(matrix_ALT)
  
  # slow
  matrix_VAF <- make_VAF_matrix(matrix_REF = matrix_REF, matrix_ALT = matrix_ALT, fast_cutoff = 0.9) # slow
  rownames(matrix_VAF) <- rownames(matrix_ALT)
  colnames(matrix_VAF) <- colnames(matrix_ALT)
  write.csv(matrix_VAF, "matrix_VAF.csv", row.names = T)
  write.csv(matrix_REF, "matrix_REF.csv", row.names = T)
  write.csv(matrix_ALT, "matrix_ALT.csv", row.names = T)
  
  # color
  anno_color <- list(celltype = celltypes_color)
  
  # plot
  matrix_VAF <- read.csv("matrix_VAF.csv", row.names = 1)
  matrix <- matrix_VAF
  matrix[is.na(matrix)] <- 0
  matrix[matrix>=2] <- 2
  graphics.off()
  p1 <- pheatmap(matrix, annotation_col = meta[,c("celltype","celltype")], cluster_cols = T, annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F,color = colorRampPalette(brewer.pal(9,"GnBu"))(50), clustering_distance_cols = "manhattan")
  p2 <- pheatmap(matrix[,meta$cell[order(meta$celltype)]], annotation_col = meta[,c("celltype","celltype")], cluster_cols = F, annotation_colors = anno_color,
                 show_rownames = F, show_colnames = F,color = colorRampPalette(brewer.pal(9,"GnBu"))(50))
  p <- plot_grid(p1$gtable, p2$gtable, ncol = 2)
  ggsave("heatmap_VAF.pdf",p, width = 20, height = 10)
}










# statistics
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
sample7 <- c("P1","P2","P3","P4","P5","P6","P8")
cutoff_perc <- c(0.9,0.9,0.9,0.9,0.9,0.9,0.9)
# color
celltypes <- paste0("type_",c(0:11))
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(8))
show_col(celltypes_color)
celltypes_color <- celltypes_color[c(1:5,13,12,9,11,6,7,10)]
names(celltypes_color) <- celltypes
for(i in c(3,7)){
  ################################################################################
  path = "H:/scclone/PMID38570491_ST"
  rna_path = paste0(path,"/paper/ST/",sample7[i])
  mut_path = paste0(path,"/",sample7[i],"/annovar")
  work_path = paste0(path,"/",sample7[i],"/work")
  setwd(paste0(work_path,"/4_vaf"))
  ################################################################################
  
  
  # read data
  single <- read.csv(paste0(work_path, "/3_filter/single_filter_more.csv"))
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  
  # plot signature
  matrix_VAF <- read.csv("matrix_VAF.csv", row.names = 1)
  single_keep <- single %>% filter(var %in% rownames(matrix_VAF))
  signature <- snv_to_sig(dat = single_keep[,c(1,2,3,4)], ref = "hg38")
  
  p <- show_signature(array = as.array(table(signature$signature)), name = "VAF")
  ggsave("signature_VAF.pdf",p$plot, width = 8, height = 3)
  
  # VAF
  matrix_ALT <- read.csv("matrix_ALT.csv", row.names = 1)
  matrix_REF <- read.csv("matrix_REF.csv", row.names = 1)
  matrix_VAF[is.na(matrix_VAF)] <- 0
  matrix_VAF[matrix_VAF>2] <- 2
  single$VAF <- single$AD_2/(single$AD_1+single$AD_2)
  matrix_VAF_raw <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","VAF")], var~cell, value.var = "VAF", fun.aggregate=sum)))
  matrix_VAF_raw[is.na(matrix_VAF_raw)] <- 0
  matrix_VAF_raw <- matrix_VAF_raw[rownames(matrix_VAF),colnames(matrix_VAF)]
  data_tem <- data.frame(VAF = as.numeric(unlist(matrix_VAF_raw)),
                         ALT = as.numeric(unlist(matrix_ALT)),
                         REF = as.numeric(unlist(matrix_REF)),
                         genotype = as.numeric(unlist(matrix_VAF)))
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
  
  p <- gg.gap(plot=p1,segments=list(c(10000,20000)),ylim=c(0,7*10^4))
  #add.legend(plot = p1,margin = c(top=1,right=1,bottom=200,left=300))
  ggsave("hist_VAF_raw.pdf",p, width = 6, height = 5)
  
  p <- gg.gap(plot=p2,segments=list(c(10000,20000)),ylim=c(0,7*10^4))
  #add.legend(plot = p2,margin = c(top=1,right=1,bottom=200,left=300))
  ggsave("hist_VAF.pdf",p, width = 6, height = 5)
  
  # plot distance
  meta_order <- meta %>% arrange(celltype)
  rownames(meta_order) <- meta_order$cell
  matrix_order <- matrix_VAF[,rownames(meta_order)]
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
  anno_color <- list(celltype = celltypes_color)
  p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("celltype","sample")], annotation_colors = anno_color,
                          annotation_row = meta_order[,c("celltype","sample")],
                          show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                          color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
  ggsave("distance.pdf", p, width = 10, height = 10)
  
  # statistic
  demo <- as.matrix(dist(t(matrix_VAF), method = "manhattan", diag = FALSE, upper = FALSE))
  hist(unlist(demo), xlim = c(0,max(demo)), breaks = 100)
  cutoff <- quantile(unlist(demo),cutoff_perc[i])
  cluster <- get_hclust(matrix = matrix_VAF, cut = cutoff, method = "manhattan")
  table(cluster)
  meta$cluster_tem <- paste0("cluster_",cluster[match(meta$cell, names(cluster))])
  #index <- get_assess_index(meta_order = meta_order, matrix_order = matrix_order)
  #write.csv(index, paste0("index_final.csv"), row.names = F)
  p <- pheatmap(matrix_VAF[,meta$cell[order(meta$cluster_tem)]], annotation_col = meta[,c("cluster_tem","celltype","sample")], cluster_cols = T, annotation_colors = anno_color,
                show_rownames = F, show_colnames = F,color = colorRampPalette(brewer.pal(9,"GnBu"))(50), clustering_distance_cols = "manhattan")
  ggsave("heatmap_VAF_cluster_tem.pdf",p, width = 10, height = 10)
}





