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
celltypes <- c("B_cell","Ductal_cell","Fibroblasts","Mac_Mono","Mast_cell","MKI67+cell","NK","Plasma_cell","T_cell")
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
      matrix_cell <- as.matrix(data@assays$RNA$counts[,idx])
      anno_cell <- data.frame(row.names = colnames(matrix_cell), anno = data$cluster_new[idx]) %>% arrange(anno)
      gene_anno <- read.table("H:/scclone/GENEpos/hg19_genes_pos.txt", sep = "\t", quote = "", header = F, row.names = 1)
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






