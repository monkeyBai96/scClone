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
  
  table <- as.data.frame(table(meta$celltype, meta$cluster_new)) %>% filter(Freq>0)
  table <- as.data.frame(table(table$Var1))
  celltype_keep <- as.character(table$Var1[which(table$Freq>1)])
  
  dir.create(paste0(work_path, "/8_figure/infercnv"), showWarnings = F)
  if(length(celltype_keep)>0){
    for(j in 1:length(celltype_keep)){
      dir.create(paste0(work_path, "/8_figure/infercnv/",celltype_keep[j]), showWarnings = F)
      setwd(paste0(work_path, "/8_figure/infercnv/",celltype_keep[j]))
      delete <- table(data$cluster_new[which(data$celltype==celltype_keep[j])])
      delete <- delete[delete>1]
      idx <- which(data$celltype==celltype_keep[j] & data$cluster_new %in% names(delete))
      matrix_cell <- as.matrix(data@assays$Spatial$counts[,idx])
      anno_cell <- data.frame(row.names = colnames(matrix_cell), anno = data$cluster_new[idx]) %>% arrange(anno)
      gene_anno <- read.table("H:/scclone/GENEpos/hg38_genes_pos.txt", sep = "\t", quote = "", header = F, row.names = 1)
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = matrix_cell,
                                           annotations_file = anno_cell,
                                           gene_order_file = gene_anno,
                                           ref_group_names = anno_cell$anno[1],
                                           delim="\t")
      options("Seurat.object.assay.version" = "v3")
      infercnv_obj = infercnv::run(infercnv_obj,
                                   cutoff=1, 
                                   out_dir=getwd(),
                                   cluster_by_groups=T, 
                                   analysis_mode="subclusters",
                                   denoise=TRUE,
                                   HMM=TRUE,
                                   num_threads=8,
                                   write_expr_matrix = T,
                                   write_phylo = T
      )
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
  
  meta$barcode <- substring(meta$file,4,nchar(meta$file))
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(meta$barcode,colnames(data))] <- meta$cluster_new
  data$cluster_new[is.na(data$cluster_new)] <- "un_labeled"
  data$celltype <- "un_labeled"
  data$celltype[match(meta$barcode,colnames(data))] <- meta$celltype
  data$sample <- sample2[i]
  
  dir.create(paste0(work_path, "/8_figure2/infercnv"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure2/infercnv"))
  matrix_cell <- as.matrix(data@assays$Spatial$counts)
  anno_cell <- data.frame(row.names = colnames(matrix_cell), anno = data$cluster_new) %>% arrange(anno)
  gene_anno <- read.table("H:/scclone/GENEpos/hg38_genes_pos.txt", sep = "\t", quote = "", header = F, row.names = 1)
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = matrix_cell,
                                       annotations_file = anno_cell,
                                       gene_order_file = gene_anno,
                                       ref_group_names = "un_labeled",
                                       delim="\t")
  options("Seurat.object.assay.version" = "v3")
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=1, 
                               out_dir=getwd(),
                               cluster_by_groups=T, 
                               analysis_mode="subclusters",
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads=8,
                               write_expr_matrix = T,
                               write_phylo = T
  )
}











