# HCC1259
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
################################################################################
path = "H:/scclone/PMID_SU_HCC"
rna_path = paste0(path,"/paper")
setwd(paste0(path,"/0_exp"))
################################################################################


#################### normalize basic scRNA-seq data ############################
count <- read.table(paste0(rna_path,"/HCC1-2-5-9_count_with_ERCC.txt"), sep = "\t", header = T)
meta <- read_excel(paste0(rna_path,"/meta_su.xlsx"))

# workflow
count <- count[!duplicated(count[,1]),]
count <- col1_to_rowname(count)
count <- count[,meta$Cell_ID]
obj <- CreateSeuratObject(counts = count, min.cells = 3)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(obj)
obj@meta.data$sample <- substring(rownames(obj@meta.data),1,4)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:10)
obj <- RunTSNE(obj, dims = 1:10)

graphics.off()
pdf("exp_umap.pdf", width = 7, height = 13)
p1 <- DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
plot_grid(p1, p2, ncol = 1)
dev.off()
pdf("exp_gene_umap.pdf", width = 11, height = 12)
FeaturePlot(object = obj, features = c("HNF4A","CD24","ANPEP","SOX9","GPC3","AFP","EPCAM","LGR5","AXIN2","COL6A1","REG3A","GPC4"), ncol = 3, reduction = "umap", pt.size = 0.1,
            cols = colorRampPalette(c("grey","royalblue","navy"))(10))
dev.off()
pdf("exp_tsne.pdf", width = 7, height = 13)
p1 <- DimPlot(object = obj, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = obj, reduction = 'tsne',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
plot_grid(p1, p2, ncol = 1)
dev.off()
pdf("exp_gene_tsne.pdf", width = 11, height = 12)
FeaturePlot(object = obj, features = c("HNF4A","CD24","ANPEP","SOX9","GPC3","AFP","EPCAM","LGR5","AXIN2","COL6A1","REG3A","GPC4"), ncol = 3, reduction = "tsne", pt.size = 0.1,
            cols = colorRampPalette(c("grey","royalblue","navy"))(10))
dev.off()

# su result
obj$seurat_clusters_su <- as.factor(meta$cluster)
embeddings <- as.matrix(data.frame(row.names = meta$Cell_ID, tsne1 = meta$tSNE1, tsne2 = meta$tSNE2))
obj[["tsne_su"]] <- CreateDimReducObject(embeddings = embeddings, key = "su_", assay = "RNA")
Idents(obj) <- "seurat_clusters_su"
celltype_su <- c("Tumor_cell","Others_2","Tumor_cell","Tumor_cell","T_cell","Tumor_cell","Macrophage","B_cell","T_cell","Others_10")
names(celltype_su) <- levels(obj)
obj <- RenameIdents(obj, celltype_su)
obj$celltype_su <- Idents(obj)
pdf("exp_tsne_su.pdf", width = 7, height = 20)
p1 <- DimPlot(object = obj, reduction = 'tsne_su',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = obj, reduction = 'tsne_su',label = TRUE, group.by = "seurat_clusters_su", pt.size = 3,raster = TRUE)
p3 <- DimPlot(object = obj, reduction = 'tsne_su',label = TRUE, group.by = "celltype_su", pt.size = 3,raster = TRUE)
plot_grid(p1, p2, p3, ncol = 1)
dev.off()
meta_exp <- obj@meta.data
meta <- data.frame(row.names = rownames(meta_exp), cell = rownames(meta_exp), barcode = rownames(meta_exp), file = rownames(meta_exp),
                   sample = meta_exp$sample, cluster_bai = meta_exp$seurat_clusters, cluster_su = meta_exp$seurat_clusters_su,
                   celltype = meta_exp$celltype_su)
write.csv(meta, "meta.csv", row.names = T)
write.csv(meta_exp, "meta_exp.csv", row.names = T)
saveRDS(obj, "HCC1259.rds")


