# ST
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
for(i in 1:2){
  ################################## path ########################################
  path = "H:/scclone/PMID32579974_ST"
  rna_path = paste0(path,"/paper/ST/GSE144239_RAW/",sample2[i])
  mut_path = paste0(path,"/",sample2[i],"/annovar")
  work_path = paste0(path,"/",sample2[i],"/work")
  setwd(paste0(work_path,"/0_exp"))
  ################################################################################
  
  matrix = Read10X(data.dir = rna_path)
  matrix = matrix[-grep("^MT-", rownames(matrix)),] #
  matrix = matrix[-grep("^RP[SL]", rownames(matrix)),] #
  data = CreateSeuratObject(counts = matrix, project = sample2[i], assay = "Spatial")
  img = Read10X_Image(paste0(rna_path,"/spatial"), image.name = "tissue_hires_image.png")
  img = img[Cells(x = data)]
  data@images$slice1 = img
  data@images$slice1@assay = "Spatial"
  data@images$slice1@key = "slice1_"
  data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data = FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data = ScaleData(data, features = rownames(data))
  data <- RunPCA(data, assay = "Spatial", verbose = FALSE)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
  data <- FindClusters(data, verbose = FALSE, resolution = 0.8)
  data <- FindClusters(data, verbose = FALSE)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30)
  
  p1 <- DimPlot(data, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 1)
  p3 <- SpatialDimPlot(data, label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.1)
  ggsave("umap_spatial.pdf",p1+p2+p3, width = 15, height = 5)
  
  data <- FindClusters(data, verbose = FALSE, resolution = 1.35)
  p1 <- DimPlot(data, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 1)
  p3 <- SpatialDimPlot(data, label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.1)
  ggsave("umap_spatial_resolution_1.5.pdf",p1+p2+p3, width = 15, height = 5)
  saveRDS(data, "st.rds")
  
  meta <- data@meta.data
  write.csv(meta, "meta_exp.csv", row.names = T)
  write.table(colnames(data@assays$Spatial$counts), paste0(sample2[i],"_barcode.txt"), col.names = F, row.names = F, quote = F)
  meta$cell <- paste0(meta$orig.ident,"_",substring(rownames(meta),1,16))
  meta$file <- paste0("CB_",rownames(meta))
  rownames(meta) <- meta$cell
  meta$celltype <- paste0("type_",meta$Spatial_snn_res.0.8)
  meta$celltype_resolution <- paste0("type_",meta$Spatial_snn_res.1.35)
  
  write.csv(meta, paste0(work_path,"/0_exp/meta.csv"), row.names = T)
  #SpatialDimPlot(data[,which(data$Spatial_snn_res.1%in%c(4,8))], label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 0.2, group.by = "Spatial_snn_res.1")
  #x <- data.frame(cell = data@images$slice1@boundaries$centroids@cells,
  #                coor = data@images$slice1@boundaries$centroids@coords,
  #                cluster = data$Spatial_snn_res.1)
  #xx <- x %>% filter(cluster %in% c(4,8))
  #xx$label <- c(1:125)
  #ggplot(xx, aes(y = -coor.x, x = coor.y, color = cluster))+
  #  geom_point()+
  #  scale_color_brewer(palette = "Set1")+
  #  geom_text_repel(aes(label = label), size = 5)
  #change <- xx %>% filter(label %in% c(39,4,108,50,51,62,111))
  #change$cell <- paste0("P4_1_",change$cell)
  #write.csv(change, "disorder_cell.csv", row.names = F)
}






