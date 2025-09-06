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
  
  table <- as.data.frame(table(meta$celltype, meta$cluster_new)) %>% filter(Freq>0)
  table <- as.data.frame(table(table$Var1))
  celltype_keep <- as.character(table$Var1[which(table$Freq>1)])
  
  dir.create(paste0(work_path, "/8_figure/subcluster"), showWarnings = F)
  if(length(celltype_keep)>0){
    for(j in 1:length(celltype_keep)){
      dir.create(paste0(work_path, "/8_figure/subcluster/",celltype_keep[j]), showWarnings = F)
      setwd(paste0(work_path, "/8_figure/subcluster/",celltype_keep[j]))
      subdata <- data[,which(colnames(data) %in% meta$barcode & data$celltype==celltype_keep[j])]
      Idents(subdata) <- "cluster_new"
      data.markers <- FindAllMarkers(subdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      saveRDS(data.markers,"data.markers.rds")
      write.csv(data.markers, "data.markers.csv", row.names = T)
      
      # by subcluster
      data.markers <- read.csv("data.markers.csv", row.names = 1)
      cluster <- as.character(unique(data.markers$cluster))
      list <- list()
      for(k in 1:length(cluster)){
        gene <- data.markers$gene[which(data.markers$cluster==cluster[k] & data.markers$avg_log2FC>0.5 & data.markers$p_val_adj<0.05)]
        if(length(gene)>5){
          gene.df <- bitr(gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
          ego <- enrichGO(gene = gene.df$ENTREZID,
                          OrgDb=org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          minGSSize = 1,
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          readable = TRUE)
          res_go <- as.data.frame(ego@result)
          res_go <- res_go[which(res_go$pvalue<0.05),]
          R.utils::setOption("clusterProfiler.download.method","auto")
          kk <- enrichKEGG(gene = gene.df$ENTREZID,
                           organism ="human",
                           keyType = "kegg",
                           pvalueCutoff = 1,
                           qvalueCutoff = 1,
                           use_internal_data = FALSE)
          res_kegg <- as.data.frame(kk@result)
          res_kegg <- res_kegg[which(res_kegg$pvalue<0.05),]
          for(m in 1:nrow(res_kegg)){
            str <- do.call(c,strsplit(res_kegg$geneID[m],split = "/"))
            str <- gene.df$SYMBOL[match(str,gene.df$ENTREZID)]
            str <- paste(str,collapse = "/")
            res_kegg$geneID[m] <- str
          }
          ench <- rbind(res_go,res_kegg)
          ench$cluster <- cluster[k]
          rownames(ench) <- NULL
          list[[k]] <- ench
        }
      }
      ench_all <- do.call(rbind,list)
      write.csv(ench_all, "enrichment.csv", row.names = F)
    }
  }
}


