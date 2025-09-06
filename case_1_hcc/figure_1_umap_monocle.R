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
path = "H:/scclone/PMID_SU_HCC"
data_all <- readRDS(paste0(path, "/0_exp/HCC1259.rds"))
for(i in 2:3){
  ################################## path ########################################
  rna_path = paste0(path,"/paper")
  mut_path = paste0(path,"/",HCC_name[i],"/annovar")
  work_path = paste0(path,"/",HCC_name[i],"/work")
  setwd(paste0(work_path, "/8_figure"))
  ################################################################################
  
  
  data <- data_all[,which(data_all$sample==HCC_name[i])]
  meta <- read.csv("meta_final.csv", row.names = 1)
  matrix <- read.csv("matrix_final.csv", row.names = 1)
  matrix <- matrix[,rownames(meta)]
  
  # color
  n = length(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell", "un_labeled"))
  color <- c(pal_npg("nrc")(10))[c(1:n)]
  show_col(color)
  names(color) <- c(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell", "un_labeled"))
  celltype_color <- color
  cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
  cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")
  
  # compare UMAP & TSNE
  umap <- umap(t(matrix), n_components = 2)
  umap_result <- as.data.frame(umap$layout) %>% mutate(celltype = meta$celltype, cluster_new = meta$cluster_new, cluster_su = paste0("seurat_",meta$cluster_su))
  colnames(umap_result)[1:2] <- c("UMAP1", "UMAP2")
  p1 <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = cluster_new))+
    geom_point()+
    scale_color_manual(values = cluster_new_color)+
    theme_base()
  p2 <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = celltype))+
    geom_point()+
    scale_color_manual(values = celltype_color)+
    theme_base()
  p3 <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = cluster_su))+
    geom_point()+
    theme_base()
  p <- plot_grid(p1, p2, p3, ncol = 1)
  ggsave("umap_var.pdf", p, width = 7, height = 16)
  
  tsne <- Rtsne(t(matrix), dims = 2)
  tsne_result <- as.data.frame(tsne$Y) %>% mutate(celltype = meta$celltype, cluster_new = meta$cluster_new, cluster_su = paste0("seurat_",meta$cluster_su))
  colnames(tsne_result)[1:2] <- c("tSNE1", "tSNE2")
  p1 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = cluster_new))+
    geom_point()+
    scale_color_manual(values = cluster_new_color)+
    theme_base()
  p2 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = celltype))+
    geom_point()+
    scale_color_manual(values = celltype_color)+
    theme_base()
  p3 <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = cluster_su))+
    geom_point()+
    theme_base()
  p <- plot_grid(p1, p2, p3, ncol = 1)
  ggsave("tsne_var.pdf", p, width = 7, height = 16)
  
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(rownames(meta),colnames(data))] <- meta$cluster_new
  p1 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "seurat_clusters_su", pt.size = 2,raster = F)
  p2 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "celltype_su", pt.size = 2,raster = F, cols = celltype_color)
  p3 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 2,raster = F)
  p4 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "cluster_new", pt.size = 2,raster = F, cols = cluster_new_color)
  p <- plot_grid(p1, p2, p3, p4, ncol = 1)
  ggsave("umap_exp.pdf", p, width = 7, height = 21)
  p1 <- DimPlot(object = data, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters_su", pt.size = 2,raster = F)
  p2 <- DimPlot(object = data, reduction = 'tsne',label = TRUE, group.by = "celltype_su", pt.size = 2,raster = F, cols = celltype_color)
  p3 <- DimPlot(object = data, reduction = 'tsne',label = TRUE, group.by = "sample", pt.size = 2,raster = F)
  p4 <- DimPlot(object = data, reduction = 'tsne',label = TRUE, group.by = "cluster_new", pt.size = 2,raster = F, cols = cluster_new_color)
  p <- plot_grid(p1, p2, p3, p4, ncol = 1)
  ggsave("tsne_exp.pdf", p, width = 7, height = 21)
  
  # monocle
  dir.create(paste0(work_path, "/8_figure/monocle"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/monocle"))
  subdata <- data[,which(colnames(data) %in% meta$cell)]
  matrix <- as(as.matrix(subdata@assays$RNA$counts),"sparseMatrix")
  pd <- new('AnnotatedDataFrame', data = subdata@meta.data)
  fData <- data.frame(gene_short_name = row.names(matrix), row.names = row.names(matrix))
  fd <- new('AnnotatedDataFrame', data = fData)
  cds <- newCellDataSet(matrix, phenoData = pd, featureData = fd,expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds  <- estimateDispersions(cds, cores=4, relative_expr = TRUE)
  cds <- detectGenes(cds, min_expr = 10) 
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.1)
  cds <- setOrderingFilter(cds, ordering_genes$gene_id)
  cds <- reduceDimension(cds, reduction_method = 'DDRTree')
  cds <- orderCells(cds, reverse = F)
  p1 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + scale_color_manual(values = cluster_new_color)
  p2 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "State")+ scale_color_brewer(palette="Set1")
  p3 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "Pseudotime")
  p4 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "celltype_su") + scale_color_manual(values = celltype_color)
  p5 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters")
  p6 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters_su")
  p <- plot_grid(p1,p2,p3,p4,p5,p6, ncol = 6)
  ggsave("all_monocle.pdf",p, width = 20, height = 4)
  n_state <- length(unique(cds$State))
  n_cluster <- length(unique(cds$cluster_new))
  p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = cluster_new_color)
  ggsave("all_monocle_split.pdf",p, width = 2.5*n_state, height = 4)
  p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~cluster_new", nrow = 1) + scale_color_manual(values = cluster_new_color)
  ggsave("all_monocle_split1.pdf",p, width = 2.5*n_cluster, height = 4)
  data_plot <- cds@phenoData@data[,c("State","cluster_new")]
  data_plot <- as.data.frame(table(data_plot)) %>% mutate(perc = NA)
  sum <- data_plot %>% group_by(State) %>% summarise(sum = sum(Freq))
  for(j in 1:nrow(data_plot)){
    data_plot$perc[j] <- data_plot$Freq[j]/sum$sum[which(sum$State==data_plot$State[j])]
  }
  p <- ggplot(data_plot, aes(x = State, y = perc, fill = cluster_new))+
    geom_bar(stat = 'identity',width = 0.9)+
    scale_fill_manual(values = cluster_new_color)+
    theme_bw(base_size = 12)
  ggsave("all_monocle_barplot.pdf",p, width = 1+0.5*n_state, height = 5)
  
  table <- as.data.frame(table(meta$celltype, meta$cluster_new)) %>% filter(Freq>0)
  table <- as.data.frame(table(table$Var1))
  celltype_keep <- as.character(table$Var1[which(table$Freq>1)])
  if(length(celltype_keep)>0){
    for(j in 1:length(celltype_keep)){
      subdata <- data[,which(colnames(data) %in% meta$cell & data$celltype_su==celltype_keep[j])]
      id <- which(rowSums(subdata@assays$RNA$counts)>ncol(subdata))
      matrix <- as(as.matrix(subdata@assays$RNA$counts[id,]),"sparseMatrix")
      pd <- new('AnnotatedDataFrame', data = subdata@meta.data)
      fData <- data.frame(gene_short_name = row.names(matrix), row.names = row.names(matrix))
      fd <- new('AnnotatedDataFrame', data = fData)
      cds <- newCellDataSet(matrix, phenoData = pd, featureData = fd,expressionFamily = negbinomial.size())
      cds <- estimateSizeFactors(cds)
      cds  <- estimateDispersions(cds, cores=4, relative_expr = TRUE)
      cds <- detectGenes(cds, min_expr = 10) 
      disp_table <- dispersionTable(cds)
      ordering_genes <- subset(disp_table, mean_expression >= 0.1)
      cds <- setOrderingFilter(cds, ordering_genes$gene_id)
      cds <- reduceDimension(cds, reduction_method = 'DDRTree')
      cds <- orderCells(cds, reverse = F)
      p1 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + scale_color_manual(values = cluster_new_color)
      p2 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "State")+ scale_color_brewer(palette="Set1")
      p3 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "Pseudotime")
      p4 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "celltype_su") + scale_color_manual(values = celltype_color)
      p5 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters")
      p6 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters_su")
      p <- plot_grid(p1,p2,p3,p4,p5,p6, ncol = 6)
      ggsave(paste0(celltype_keep[j],"_monocle.pdf"),p, width = 20, height = 4)
      n_state <- length(unique(cds$State))
      n_cluster <- length(unique(cds$cluster_new))
      p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = cluster_new_color)
      ggsave(paste0(celltype_keep[j],"_monocle_split.pdf"),p, width = 2.5*n_state, height = 4)
      p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~cluster_new", nrow = 1) + scale_color_manual(values = cluster_new_color)
      ggsave(paste0(celltype_keep[j],"_monocle_split1.pdf"),p, width = 2.5*n_cluster, height = 4)
      data_plot <- cds@phenoData@data[,c("State","cluster_new")]
      data_plot <- as.data.frame(table(data_plot)) %>% mutate(perc = NA)
      sum <- data_plot %>% group_by(State) %>% summarise(sum = sum(Freq))
      for(k in 1:nrow(data_plot)){
        data_plot$perc[k] <- data_plot$Freq[k]/sum$sum[which(sum$State==data_plot$State[k])]
      }
      p <- ggplot(data_plot, aes(x = State, y = perc, fill = cluster_new))+
        geom_bar(stat = 'identity',width = 0.9)+
        scale_fill_manual(values = cluster_new_color)+
        theme_bw(base_size = 12)
      ggsave(paste0(celltype_keep[j],"_monocle_barplot.pdf"),p, width = 1+0.5*n_state, height = 5)
    }
  }
}



