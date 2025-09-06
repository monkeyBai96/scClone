# PMID31558476
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
################################## path ########################################
path = "H:/scclone/PMID31558476_bm"
rna_path = paste0(path,"/paper")
work_path = paste0(path,"/work")
setwd(paste0(work_path,"/0_exp"))
################################################################################


#################### normalize basic scRNA-seq data ############################
# read meta
meta1 <- read.table(paste0(rna_path,"/GSE110499_SraRunTable.txt"), sep = ",", header = T) %>% filter(GEO_Accession..exp. != "GSM2994860")
meta2 <- read.table(paste0(rna_path,"/GSE106218_SraRunTable.txt"), sep = ",", header = T) %>% filter(facs_sorting != "CD138 negative")
meta <- rbind(meta1[,c("Run","Sample.Name","prep.site","source_name")], meta2[,c("Run","Sample.Name","prep.site","source_name")])
colnames(meta) <- c("SRR","Sample_id","tissue","sample")
cell_name <- read.table(paste0(rna_path, "/cell_name.txt"), sep = "\t")
meta$cell <- cell_name$V2[match(meta$Sample_id, cell_name$V1)]
meta <- meta %>% filter(sample %in% c("MM02","MM02EM","MM17","MM25","MM30","MM34","MM34EM","MM36EM","MM16","MM16R"))
# read data
raw1 <- read.table(paste0(rna_path,"/GSE110499_GEO_processed_MM_raw_TPM_matrix.txt"), row.names = 1, header = T)
raw2 <- read.table(paste0(rna_path,"/GSE106218_GEO_processed_MM_raw_TPM_matrix.txt"), row.names = 1, header = T)
raw <- cbind(raw1, raw2)
raw <- raw[,!duplicated(colnames(raw))]
colnames(raw) <- gsub("_0","_",colnames(raw))
meta$cell <- gsub("_0","_",meta$cell)
setdiff(colnames(raw), meta$cell)
setdiff(meta$cell, colnames(raw))

raw <- raw[,which(colnames(raw) %in% meta$cell)]
meta <- meta %>% filter(cell %in% colnames(raw))
rownames(meta) <- meta$cell

## scRNA workflow
obj <- CreateSeuratObject(counts = raw, min.cells = 5, meta.data = data.frame(row.names = meta$cell, sample = meta$sample, tissue = meta$tissue))
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(obj)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, dims = 1:10)
obj <- RunTSNE(obj, dims = 1:10)
DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)

# save
graphics.off()
p1 <- DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
p3 <- DimPlot(object = obj, reduction = 'umap',label = TRUE, group.by = "tissue", pt.size = 3,raster = TRUE)
p <- plot_grid(p1, p2, p3, ncol = 1)
ggsave("exp_umap.pdf",p, width = 7, height = 19)
p1 <- DimPlot(object = obj, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters", pt.size = 3,raster = TRUE)
p2 <- DimPlot(object = obj, reduction = 'tsne',label = TRUE, group.by = "sample", pt.size = 3,raster = TRUE)
p3 <- DimPlot(object = obj, reduction = 'tsne',label = TRUE, group.by = "tissue", pt.size = 3,raster = TRUE)
p <- plot_grid(p1, p2, p3, ncol = 1)
ggsave("exp_tsne.pdf",p, width = 7, height = 19)
saveRDS(obj, "PMID31558476.rds")
write.csv(meta, "meta_exp.csv", row.names = T)
meta <- meta %>% mutate(file = SRR, barcode = cell)
meta$celltype <- paste0("type_", obj$seurat_clusters[match(meta$cell, colnames(obj))])
write.csv(meta, "meta.csv", row.names = T)




