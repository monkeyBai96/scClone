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
sample2 <- c("P3_YF","P3_ZY")
################################################################################
path = "H:/scclone/PMID37612267"
rna_path = paste0(path,"/paper/single")
setwd(paste0(path,"/0_exp"))
################################################################################


#################### normalize basic scRNA-seq data ############################
for(i in 1:2){
  # import data
  data <- Read10X(data.dir = paste0(rna_path,"/",sample2[i]))
  data <- CreateSeuratObject(counts = data, project = sample2[i], min.cells = 3, min.features = 200)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") # Calculate percent mitochondrial genes
  data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[sl]") # Calculate percent ribosome genes
  data[["percent.hb"]] <- PercentageFeatureSet(data, pattern = "^HB[a-b]-") # Calculate the percent hemoglobin genes
  data <- subset(data, subset = percent.mt < 20 & nFeature_RNA > 500 & nFeature_RNA < 5000)
  # QC
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  data <- ScaleData(data, features = rownames(data), verbose = FALSE)
  data$sample <- sample2[i]
  DefaultAssay(data) <- 'RNA'  
  data <- RunPCA(data, features = VariableFeatures(object = data), verbose = FALSE)
  saveRDS(data, paste0(rna_path,"/",sample2[i],".rds"))
}
data1 <- readRDS(paste0(rna_path,"/",sample2[1],".rds"))
data2 <- readRDS(paste0(rna_path,"/",sample2[2],".rds"))

features <- SelectIntegrationFeatures(object.list = list(data1,data2), nfeatures = 10000)
anchors <- FindIntegrationAnchors(object.list = list(data1,data2), dims = 1:30, anchor.features = features)
data.combined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(data.combined) <- "integrated"
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5) # 27 cluster
data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:30)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(data.markers,paste0(path,"/0_exp/data.markers.rds"))


markers <- c("KLRD1","NKG7","GZMH","SLC30A1","CCL5","CD3E","CD3D","GZMB","GZMA","TRAC","IL2RA","IL7R","JUN",
             "FOS","CENPF","MKI67","ASPM","TOP2A","MMP7","TSPAN8","TFF2","AGR2","KRT19","CEACAM5","MUC1","EPCAM",
             "FXYD3","CEACAM6","LCN2","FXYD2","ELF3","CELA3A","CELA3B","CPA1","LYZ","CELA2A","SST","INS","PRSS1",
             "COL1A2","DCN","COL1A1","LUM","ACTA2","COL3A1","PECAM1","PLVAP","FCER1G","CD68","APOE","CXCL2","C1QA",
             "CXCL3","HLA−DRA","CD74","AIF1","CXCL1","CCL2","S100A9","S100A8","CLU","TPSB2","MS4A2","CD69","CD79A",
             "MS4A1","TNFRSF13C","LY9","BANK1","IGHA1","SSR4","MZB1","FKBP")
pdf(paste0(path,"/0_exp/marker_dot.pdf"), width = 15, height = 8)
DotPlot(data.combined, features = markers) + RotatedAxis()
dev.off()

celltype <- data.frame(id = c(0:15),
                       name = c("NK","NK","Ductal_cell","Mac/Mono","Mac/Mono","T_cell","Fibroblasts","NK","NK","B_cell",
                                "Mast_cell","Ductal_cell","Fibroblasts","Ductal_cell","MKI67+cell","Plasma_cell"))
data.combined$celltype <- celltype$name[match(data.combined$seurat_clusters,celltype$id)]

# change path
graphics.off()
pdf("exp_umap.pdf", width = 7, height = 18)
p1 <- DimPlot(object = data.combined, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = data.combined, reduction = 'umap',label = TRUE, group.by = "celltype", pt.size = 3,raster = TRUE)
p3 <- DimPlot(object = data.combined, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
plot_grid(p1, p2, p3, ncol = 1)
dev.off()
pdf("exp_tsne.pdf", width = 7, height = 18)
p1 <- DimPlot(object = data.combined, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = data.combined, reduction = 'tsne',label = TRUE, group.by = "celltype", pt.size = 3,raster = TRUE)
p3 <- DimPlot(object = data.combined, reduction = 'tsne',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
plot_grid(p1, p2, p3, ncol = 1)
dev.off()
saveRDS(data.combined,paste0(path,"/0_exp/data.combined.rds"))

meta <- data.combined@meta.data
meta$cell <- paste0(meta$sample,"_",substring(rownames(meta),1,16))
meta$file <- paste0("CB_",substring(rownames(meta),1,18))
meta$barcode <- substring(meta$file,4,nchar(meta$file))
write.table(substring(rownames(meta)[which(meta$sample=="P3_YF")],1,18), paste0(path,"/0_exp/P3_YF_barcode.txt"), col.names = F, row.names = F, quote = F)
write.table(substring(rownames(meta)[which(meta$sample=="P3_ZY")],1,18), paste0(path,"/0_exp/P3_ZY_barcode.txt"), col.names = F, row.names = F, quote = F)
write.csv(meta, paste0(path,"/0_exp/meta_exp.csv"), row.names = T)
rownames(meta) <- meta$cell
write.csv(meta[,c(1:6,8,9,11:13,10,7)], "meta.csv", row.names = T)



