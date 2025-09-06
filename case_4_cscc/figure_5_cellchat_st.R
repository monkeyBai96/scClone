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
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(meta$barcode,colnames(data))] <- meta$cluster_new
  data$celltype <- "un_labeled"
  data$celltype[match(meta$barcode,colnames(data))] <- meta$celltype
  data$sample <- sample2[i]
  
  celltype <- unique(meta$celltype)
  cluster_new <- unique(meta$cluster_new)
  
  # cellchat
  dir.create(paste0(work_path, "/8_figure/cellchat"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/cellchat"))
  idx <- which(data$cluster_new!="un_labeled")
  mat <- data@assays$Spatial$counts[,idx]
  cellchat <- createCellChat(object = mat, meta = data@meta.data[idx,], group.by = "cluster_new")
  cellchat <- setIdent(cellchat, ident.use = "cluster_new")
  cellchat@idents <- factor(cellchat@idents, levels = cluster_new)
  groupSize <- as.numeric(table(cellchat@idents))
  save <- as.data.frame(CellChatDB.human$interaction)
  write.csv(save, "cellchatDB_interaction.csv", row.names = F)
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multicore", workers = 8) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 0)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, "cellchat.rds")
  
  c1 <- readRDS("cellchat.rds")
  p1 <- netVisual_circle(c1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  p2 <- netVisual_circle(c1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  p <- plot_grid(p1, p2)
  ggsave("cellchat_all.pdf",p, width = 12, height = 6)
  mat <- c1@net$weight
  n_cluster <- ncol(mat)
  plot_list <- list()
  for(j in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    plot_list[[j]] <- netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  p <- plot_grid(plotlist = plot_list, ncol = 5)
  ggsave("cellchat_all_split.pdf",p, width = 15, height = 15)
  
  dir.create(paste0(work_path, "/8_figure/cellchat/pathway"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/cellchat/pathway"))
  graphics.off()
  ways <- c1@netP$pathways
  if(!is.null(ways)){
    for(j in 1:length(ways)){
      p1 <- netVisual_aggregate(c1, signaling = ways[j], layout = "chord")
      p2 <- netVisual_bubble(c1, signaling = ways[j])
      p <- plot_grid(p1, p2, ncol = 2)
      ggsave(paste0(ways[j],".pdf"),p, width = 6, height = 3)
    }
  }
}





# spatial cellchat
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
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(meta$barcode,colnames(data))] <- meta$cluster_new
  data$celltype <- "un_labeled"
  data$celltype[match(meta$barcode,colnames(data))] <- meta$celltype
  data$sample <- sample2[i]
  
  celltype <- unique(meta$celltype)
  cluster_new <- unique(meta$cluster_new)
  
  # cellchat
  dir.create(paste0(work_path, "/8_figure/cellchat_ST"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/cellchat_ST"))
  data_input <- Seurat::GetAssayData(data, slot = "data", assay = "Spatial") 
  meta_input <- data@meta.data[,c("cluster_new","celltype")]
  spatial.locs = Seurat::GetTissueCoordinates(data, scale = NULL,cols = c("imagerow", "imagecol"))[,1:2]
  scale.factors=data@images$slice1@scale.factors
  spot.size = 65
  conversion.factor=spot.size/scale.factors$spot
  spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
  
  #idx <- which(data$cluster_new!="un_labeled")
  #cellchat <- createCellChat(object = data_input[,idx], meta = meta_input[idx,], group.by = "cluster_new",datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
  cellchat <- createCellChat(object = data_input, meta = meta_input, group.by = "cluster_new",datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multicore", workers = 8) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, scale.distance = 0.01, contact.knn.k = 6)
  cellchat <- filterCommunication(cellchat, min.cells = 2)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, "cellchat.rds")
  
  c1 <- readRDS("cellchat.rds")
  p1 <- netVisual_circle(c1@net$count, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  p2 <- netVisual_circle(c1@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  p <- plot_grid(p1, p2)
  ggsave("cellchat_all.pdf",p, width = 12, height = 6)
  mat <- c1@net$weight
  n_cluster <- ncol(mat)
  plot_list <- list()
  for(j in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    plot_list[[j]] <- netVisual_circle(mat2, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  p <- plot_grid(plotlist = plot_list, ncol = 5)
  ggsave("cellchat_all_split.pdf",p, width = 15, height = 15)
  pdf("cellchat_heatmap.pdf",width = 8, height = 8)
  netVisual_heatmap(cellchat,measure="count",color.heatmap="GnBu",font.size=18,font.size.title=28)
  dev.off()
  c1 <- netAnalysis_computeCentrality(c1,slot.name="netP")
  
  dir.create(paste0(work_path, "/8_figure/cellchat_ST/pathway"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/cellchat_ST/pathway"))
  graphics.off()
  ways <- c1@netP$pathways
  if(!is.null(ways)){
    for(j in 1:length(ways)){
      p1 <- netVisual_aggregate(c1, signaling = ways[j], layout = "chord")
      p2 <- netVisual_aggregate(c1, signaling = ways[j], layout="spatial", edge.width.max=2, vertex.size.max=1, alpha.image=0.2, vertex.label.cex=3.5,vertex.weight="incoming")
      p3 <- netVisual_bubble(c1, signaling = ways[j])
      p <- plot_grid(p1,p2,p3, ncol = 3)
      ggsave(paste0(ways[j],".pdf"),p, width = 20, height = 6)
    }
  }
}
















