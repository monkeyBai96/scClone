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
  setwd(paste0(work_path,"/8_figure"))
  #################################################################################
  
  
  # read data
  meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
  matrix_VAF <- read.csv(paste0(work_path,"/7_rf/matrix_rf.csv"), row.names = 1)
  matrix_VAF[is.na(matrix_VAF)] <- 0
  matrix_VAF[matrix_VAF>2] <- 2
  
  meta <- get_hclust_celltype(matrix = matrix_VAF, meta = meta)
  table(meta$cluster_new ,meta$celltype)
  table <- as.data.frame(table(paste0(meta$cluster_new, "_",meta$celltype)))
  meta_select <- meta %>% filter(paste0(meta$cluster_new, "_",meta$celltype) %in% table$Var1[which(table$Freq>5)]) %>% arrange(celltype, cluster_new)
  meta_select <- refine_meta_final(meta = meta_select)
  
  # color
  celltype_color <- celltypes_color
  cluster_new_color <- get_color_cluster_new(meta = meta_select, celltype_color = celltype_color)
  heatmap_color <- list(celltype = celltype_color, cluster_new = cluster_new_color)
  sankey_color <- c(celltype_color,cluster_new_color)
  
  p1 <- pheatmap(matrix_VAF[,meta_select$cell], annotation_col = meta_select[,c("cluster_new","celltype")], show_rownames = F, show_colnames = F, 
                 cluster_cols = F, cluster_rows = T, color = colorRampPalette(brewer.pal(9,"GnBu"))(50),
                 clustering_distance_cols = "manhattan", annotation_colors = heatmap_color)
  
  df <- meta_select %>% make_long(sample, celltype, cluster_new) 
  df$node <- factor(df$node, levels = c(unique(meta_select$sample),
                                        unique(meta_select$celltype),
                                        unique(meta_select$cluster_new)))
  p2 <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.7, node.color = "gray70", smooth = 5, width = 0.2) +
    geom_sankey_label(aes(fill = factor(node)), size = 5, color = "white") +
    scale_fill_manual(values = sankey_color)+
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
    ggtitle("sample-celltype-cluster")
  p2
  
  table <- as.data.frame(table(meta_select$cluster_new, meta_select$celltype))
  p3 <- ggplot(table, aes(area = Freq, fill = Var1, label = Var1, subgroup = Var2)) +
    geom_treemap() +
    geom_treemap_text(fontface = "italic", colour = "white", place = "centre",grow = TRUE,alpha=.6)+
    geom_treemap_subgroup_border(color='white')+
    scale_fill_manual(values = sankey_color)
  p3
  
  # save
  graphics.off()
  p <- plot_grid(p1$gtable, p2, p3, ncol = 3)
  ggsave("demo.pdf",p, width = 30, height = 10)
  write.csv(meta_select, "meta_final.csv", row.names = T)
  write.csv(matrix_VAF[,meta_select$cell], "matrix_final.csv")
}

