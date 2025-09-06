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
  setwd(paste0(work_path, "/8_figure"))
  ################################################################################
  
  
  # read data
  data <- readRDS(paste0(work_path,"/0_exp/st.rds"))
  meta <- read.csv("meta_final.csv",row.names = 1)
  matrix <- read.csv("matrix_final.csv",row.names = 1)
  celltype_color <- c(celltypes_color, "un_labeled"="gray")
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
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(meta$barcode,colnames(data))] <- meta$cluster_new
  data$celltype <- "un_labeled"
  data$celltype[match(meta$barcode,colnames(data))] <- meta$celltype
  data$sample <- sample2[i]
  p1 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 2,raster = F)
  p2 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "celltype", pt.size = 2,raster = F, cols = celltype_color)
  p3 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 2,raster = F)
  p4 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "cluster_new", pt.size = 2,raster = F, cols = cluster_new_color)
  p <- plot_grid(p1, p2, p3, p4, ncol = 1)
  ggsave("umap_exp.pdf", p, width = 7, height = 21)
  
  # monocle
  dir.create(paste0(work_path, "/8_figure/monocle"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/monocle"))
  subdata <- data[,which(colnames(data) %in% meta$barcode)]
  matrix <- as(as.matrix(subdata@assays$Spatial$counts),"sparseMatrix")
  pd <- new('AnnotatedDataFrame', data = subdata@meta.data)
  fData <- data.frame(gene_short_name = row.names(matrix), row.names = row.names(matrix))
  fd <- new('AnnotatedDataFrame', data = fData)
  cds <- newCellDataSet(matrix, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds  <- estimateDispersions(cds)
  cds <- detectGenes(cds, min_expr = 10) 
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.1)
  cds <- setOrderingFilter(cds, ordering_genes$gene_id)
  cds <- reduceDimension(cds, reduction_method = 'DDRTree')
  cds <- orderCells(cds, reverse = F)
  p1 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + scale_color_manual(values = cluster_new_color)
  p2 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "State")+ scale_color_brewer(palette="Set1")
  p3 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "Pseudotime")
  p4 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "celltype") + scale_color_manual(values = celltype_color)
  p5 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters")
  p <- plot_grid(p1,p2,p3,p4,p5,ncol = 5)
  ggsave("all_monocle.pdf",p, width = 17, height = 4)
  n_state <- length(unique(cds$State))
  n_cluster <- length(unique(cds$cluster_new))
  p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = cluster_new_color)
  ggsave("all_monocle_split.pdf",p, width = 2.5*n_state, height = 4,limitsize = FALSE)
  p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~cluster_new", nrow = 1) + scale_color_manual(values = cluster_new_color)
  ggsave("all_monocle_split1.pdf",p, width = 2.5*n_cluster, height = 4,limitsize = FALSE)
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
  ggsave("all_monocle_barplot.pdf",p, width = 2+1*n_state, height = 5)
  
  table <- as.data.frame(table(meta$celltype, meta$cluster_new)) %>% filter(Freq>0)
  table <- as.data.frame(table(table$Var1))
  celltype_keep <- as.character(table$Var1[which(table$Freq>1)])
  if(length(celltype_keep)>0){
    for(j in 1:length(celltype_keep)){
      subdata <- data[,which(colnames(data) %in% meta$barcode & data$celltype==celltype_keep[j])]
      matrix <- as(as.matrix(subdata@assays$Spatial$counts),"sparseMatrix")
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
      p4 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "celltype") + scale_color_manual(values = celltype_color)
      p5 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters")
      p <- plot_grid(p1,p2,p3,p4,p5,ncol = 5)
      ggsave(paste0(celltype_keep[j],"_monocle.pdf"),p, width = 17, height = 4)
      n_state <- length(unique(cds$State))
      n_cluster <- length(unique(cds$cluster_new))
      p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = cluster_new_color)
      ggsave(paste0(celltype_keep[j],"_monocle_split.pdf"),p, width = 2.5*n_state, height = 4,limitsize = FALSE)
      p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~cluster_new", nrow = 1) + scale_color_manual(values = cluster_new_color)
      ggsave(paste0(celltype_keep[j],"_monocle_split1.pdf"),p, width = 2.5*n_cluster, height = 4,limitsize = FALSE)
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
      ggsave(paste0(celltype_keep[j],"_monocle_barplot.pdf"),p, width = 2+1*n_state, height = 5)
    }
  }
}









# 2 cluster
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
  setwd(paste0(work_path, "/8_figure2"))
  ################################################################################
  
  
  # read data
  data <- readRDS(paste0(work_path,"/0_exp/st.rds"))
  meta <- read.csv("meta_final.csv",row.names = 1)
  matrix <- read.csv("matrix_final.csv",row.names = 1)
  celltype_color <- c(celltypes_color, "un_labeled"="gray")
  if(i==2){
    cluster_new_color <- c("cluster_1"=pal_nejm("default")(8)[7], "cluster_2"=pal_nejm("default")(8)[5])
  }else{
    cluster_new_color <- c("cluster_1"=pal_nejm("default")(8)[5], "cluster_2"=pal_nejm("default")(8)[7])
  }
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
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(meta$barcode,colnames(data))] <- meta$cluster_new
  data$celltype <- "un_labeled"
  data$celltype[match(meta$barcode,colnames(data))] <- meta$celltype
  data$sample <- sample2[i]
  p1 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 2,raster = F)
  p2 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "celltype", pt.size = 2,raster = F, cols = celltype_color)
  p3 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 2,raster = F)
  p4 <- DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "cluster_new", pt.size = 2,raster = F, cols = cluster_new_color)
  p <- plot_grid(p1, p2, p3, p4, ncol = 1)
  ggsave("umap_exp.pdf", p, width = 7, height = 21)
  
  # monocle
  dir.create(paste0(work_path, "/8_figure2/monocle"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure2/monocle"))
  subdata <- data[,which(colnames(data) %in% meta$barcode)]
  matrix <- as(as.matrix(subdata@assays$Spatial$counts),"sparseMatrix")
  pd <- new('AnnotatedDataFrame', data = subdata@meta.data)
  fData <- data.frame(gene_short_name = row.names(matrix), row.names = row.names(matrix))
  fd <- new('AnnotatedDataFrame', data = fData)
  cds <- newCellDataSet(matrix, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds  <- estimateDispersions(cds)
  cds <- detectGenes(cds, min_expr = 10) 
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.1)
  cds <- setOrderingFilter(cds, ordering_genes$gene_id)
  cds <- reduceDimension(cds, reduction_method = 'DDRTree')
  cds <- orderCells(cds, reverse = F)
  p1 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + scale_color_manual(values = cluster_new_color)
  p2 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "State")+ scale_color_brewer(palette="Set1")
  p3 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "Pseudotime")
  p4 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "celltype") + scale_color_manual(values = celltype_color)
  p5 <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "seurat_clusters")
  p <- plot_grid(p1,p2,p3,p4,p5,ncol = 5)
  ggsave("all_monocle.pdf",p, width = 17, height = 4)
  n_state <- length(unique(cds$State))
  n_cluster <- length(unique(cds$cluster_new))
  p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = cluster_new_color)
  ggsave("all_monocle_split.pdf",p, width = 2.5*n_state, height = 4,limitsize = FALSE)
  p <- plot_cell_trajectory(cds, markers_linear = TRUE, color_by = "cluster_new") + facet_wrap("~cluster_new", nrow = 1) + scale_color_manual(values = cluster_new_color)
  ggsave("all_monocle_split1.pdf",p, width = 2.5*n_cluster, height = 4,limitsize = FALSE)
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
  ggsave("all_monocle_barplot.pdf",p, width = 2+1*n_state, height = 5)
  
}




