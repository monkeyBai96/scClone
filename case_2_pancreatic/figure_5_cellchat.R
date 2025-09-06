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
path = "H:/scclone/PMID37612267"
data_all <- readRDS(paste0(path, "/0_exp/data.combined.rds"))
# color
celltypes <- c("B_cell","Ductal_cell","Fibroblasts","Mac/Mono","Mast_cell","MKI67+cell","NK","Plasma_cell","T_cell")
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(1:9)]
names(celltypes_color) <- celltypes
celltype_color <- c(celltypes_color, "un_labeled"="gray")
for(i in 1:1){
  ################################## path ########################################
  path = "H:/scclone/PMID37612267"
  rna_path = paste0(path,"/paper")
  mut_path = paste0(path,"/",sample2[i],"/annovar")
  work_path = paste0(path,"/",sample2[i],"/work")
  setwd(paste0(work_path,"/8_figure"))
  #################################################################################
  
  
  data <- data_all[,which(data_all$sample==sample2[i])]
  meta <- read.csv("meta_final.csv", row.names = 1)
  matrix <- read.csv("matrix_final.csv", row.names = 1)
  matrix <- matrix[,rownames(meta)]
  meta$barcode <- paste0(meta$barcode,"_",i)
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(meta$barcode,colnames(data))] <- meta$cluster_new
  meta$celltype[which(meta$celltype=="Mac/Mono")] <- "Mac_Mono"
  data$celltype[which(data$celltype=="Mac/Mono")] <- "Mac_Mono"
  
  celltype <- unique(meta$celltype)
  cluster_new <- unique(meta$cluster_new)
  
  # cellchat
  dir.create(paste0(work_path, "/8_figure/cellchat"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure/cellchat"))
  idx <- which(data$cluster_new!="un_labeled")
  mat <- data@assays$RNA$counts[,idx]
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















