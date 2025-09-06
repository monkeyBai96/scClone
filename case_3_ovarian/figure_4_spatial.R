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
  setwd(paste0(work_path, "/8_figure"))
  ################################################################################
  
  
  # read data
  data <- readRDS(paste0(work_path,"/0_exp/",sample7[i],".rds"))
  meta <- read.csv("meta_final.csv",row.names = 1)
  matrix <- read.csv("matrix_final.csv",row.names = 1)
  celltype_color <- c(celltypes_color, "un_labeled"="gray")
  cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
  cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")
  heatmap_color <- list(celltype = celltype_color, cluster_new = cluster_new_color)
  sankey_color <- c(celltype_color,cluster_new_color)
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new <- meta$cluster_new[match(colnames(data), substring(meta$file,4,21))]
  data$celltype <- "un_labeled"
  data$celltype <- meta$celltype[match(colnames(data), substring(meta$file,4,21))]
  data$sample <- sample7[i]
  
  
  celltype <- unique(meta$celltype)
  cluster_new <- unique(meta$cluster_new)
  p1 <- SpatialDimPlot(data[,!is.na(data$celltype)], label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.1, crop = F, 
                       group.by = "celltype", cols = celltype_color)
  list <- list()
  for(j in 1:length(celltype)){
    list[[j]] <- SpatialDimPlot(data[,which(data$celltype==celltype[j])], label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.1, crop = F,
                                group.by = "celltype", cols = celltype_color)
  }
  p2 <- cowplot::plot_grid(plotlist = list, ncol = length(celltype))
  ggsave(paste0("ST_celltype_all.pdf"),p1, width = 6, height = 5, limitsize = FALSE)
  ggsave(paste0("ST_celltype_split.pdf"),p2, width = 6*length(celltype), height = 5, limitsize = FALSE)
  
  p1 <- SpatialDimPlot(data[,!is.na(data$cluster_new)], label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.1, crop = F, 
                       group.by = "cluster_new", cols = cluster_new_color)
  list <- list()
  for(j in 1:length(celltype)){
    list[[j]] <- SpatialDimPlot(data[,which(data$celltype==celltype[j])], label = F, label.size = 3, image.scale = "hires", image.alpha = 0.1, crop = F,
                                group.by = "cluster_new", cols = cluster_new_color)
  }
  p2 <- cowplot::plot_grid(plotlist = list, ncol = length(celltype))
  ggsave(paste0("ST_cluster_new_all.pdf"),p1, width = 6, height = 5, limitsize = FALSE)
  ggsave(paste0("ST_cluster_new_split.pdf"),p2, width = 6*length(celltype), height = 5, limitsize = FALSE)
}














# figure2
# split 2 clone
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
  setwd(paste0(work_path,"/8_figure2"))
  ################################################################################
  
  
  # read data
  data <- readRDS(paste0(work_path,"/0_exp/",sample7[i],".rds"))
  meta <- read.csv("meta_final.csv",row.names = 1)
  matrix <- read.csv("matrix_final.csv",row.names = 1)
  celltype_color <- c(celltypes_color, "un_labeled"="gray")
  cluster_new_color <- c("cluster_1"=pal_nejm("default")(8)[5], "cluster_2"=pal_nejm("default")(8)[7])
  cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")
  heatmap_color <- list(celltype = celltype_color, cluster_new = cluster_new_color)
  sankey_color <- c(celltype_color,cluster_new_color)
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new <- meta$cluster_new[match(colnames(data), substring(meta$file,4,21))]
  data$celltype <- "un_labeled"
  data$celltype <- meta$celltype[match(colnames(data), substring(meta$file,4,21))]
  data$sample <- sample7[i]
  
  cluster_new <- unique(meta$cluster_new)
  celltype <- unique(meta$celltype)
  
  p1 <- SpatialDimPlot(data[,!is.na(data$cluster_new)], label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.5, crop = F, 
                       group.by = "cluster_new", cols = cluster_new_color)
  list <- list()
  for(j in 1:length(celltype)){
    list[[j]] <- SpatialDimPlot(data[,which(data$celltype==celltype[j])], label = F, label.size = 3, image.scale = "hires", image.alpha = 0.1, crop = F,
                                group.by = "cluster_new", cols = cluster_new_color)
  }
  p2 <- cowplot::plot_grid(plotlist = list, ncol = length(celltype))
  ggsave(paste0("ST_cluster_new_all.pdf"),p1, width = 6, height = 5)
  ggsave(paste0("ST_cluster_new_split.pdf"),p2, width = 6*length(celltype), height = 5,limitsize = FALSE)
}


















