# PMID37612267
# exp scRNA-seq workflow
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
for(i in 1:7){
  ################################################################################
  path = "H:/scclone/PMID38570491_ST"
  rna_path = paste0(path,"/paper/ST/",sample7[i])
  mut_path = paste0(path,"/",sample7[i],"/annovar")
  work_path = paste0(path,"/",sample7[i],"/work")
  setwd(paste0(work_path,"/0_exp"))
  ################################################################################
  
  
  matrix = Read10X(data.dir = rna_path)
  matrix = matrix[-grep("^MT-", rownames(matrix)),] #
  matrix = matrix[-grep("^RP[SL]", rownames(matrix)),] #
  data = CreateSeuratObject(counts = matrix, project = sample7[i], assay = "Spatial")
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
  data <- FindClusters(data, verbose = FALSE, resolution = 1)
  data <- FindClusters(data, verbose = FALSE)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30)
  
  p1 <- DimPlot(data, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3, image.scale = "hires", image.alpha = 1)
  p <- p1 + p2
  ggsave(paste0(sample7[i],"_umap_spatial.pdf"),p, width = 10, height = 5)
  saveRDS(data, paste0(sample7[i],".rds"))
  
  meta <- data@meta.data
  write.csv(meta, "meta_exp.csv", row.names = T)
  write.table(colnames(data@assays$Spatial$counts), paste0(sample7[i],"_barcode.txt"),col.names = F, row.names = F, quote = F)
  meta$cell <- paste0(meta$orig.ident,"_",substring(rownames(meta),1,16))
  meta$file <- paste0("CB_",rownames(meta))
  rownames(meta) <- meta$cell
  meta$celltype <- paste0("type_",meta$seurat_clusters)
  meta$sample <- meta$orig.ident
  write.csv(meta, "meta.csv", row.names = T)
}


















