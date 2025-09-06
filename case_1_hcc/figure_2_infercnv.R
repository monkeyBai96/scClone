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
HCC_name <- c("HCC1","HCC2","HCC5","HCC9")
path = "H:/scclone/PMID_SU_HCC"
data_all <- readRDS(paste0(path, "/0_exp/HCC1259.rds"))
for(i in 2:3){
  ################################## path ########################################
  rna_path = paste0(path,"/paper")
  mut_path = paste0(path,"/",HCC_name[i],"/annovar")
  work_path = paste0(path,"/",HCC_name[i],"/work")
  setwd(paste0(work_path, "/8_figure"))
  ################################################################################
  
  
  data <- data_all[,which(data_all$sample==HCC_name[i])]
  meta <- read.csv("meta_final.csv", row.names = 1)
  matrix <- read.csv("matrix_final.csv", row.names = 1)
  matrix <- matrix[,rownames(meta)]
  data$cluster_new <- "un_labeled"
  data$cluster_new[match(rownames(meta),colnames(data))] <- meta$cluster_new
  
  table <- as.data.frame(table(meta$celltype, meta$cluster_new)) %>% filter(Freq>0)
  table <- as.data.frame(table(table$Var1))
  celltype_keep <- as.character(table$Var1[which(table$Freq>1)])
  
  dir.create(paste0(work_path, "/8_figure/infercnv"), showWarnings = F)
  if(length(celltype_keep)>0){
    for(j in 1:length(celltype_keep)){
      dir.create(paste0(work_path, "/8_figure/infercnv/",celltype_keep[j]), showWarnings = F)
      setwd(paste0(work_path, "/8_figure/infercnv/",celltype_keep[j]))
      delete <- table(data$cluster_new[which(data$celltype_su==celltype_keep[j])])
      delete <- delete[delete>1]
      idx <- which(data$celltype_su==celltype_keep[j] & data$cluster_new %in% names(delete))
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






