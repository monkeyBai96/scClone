# 10000 vcf.g# 10000 vcf.gz
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
  data$cluster_new <- meta$cluster_new[match(colnames(data), substring(meta$file,4,21))]
  data$celltype <- "un_labeled"
  data$celltype <- meta$celltype[match(colnames(data), substring(meta$file,4,21))]
  data$sample <- sample2[i]
  
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
        if(i==1&j==1&k==2){
          gene <- data.markers$gene[which(data.markers$cluster==cluster[k] & data.markers$avg_log2FC>0.2 & data.markers$p_val<0.05)]
        }else{
          gene <- data.markers$gene[which(data.markers$cluster==cluster[k] & data.markers$avg_log2FC>0.5 & data.markers$p_val_adj<0.05)]
        }
        
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
          if(is.null(kk)==F){
            res_kegg <- as.data.frame(kk@result)
            res_kegg <- res_kegg[which(res_kegg$pvalue<0.05),]
            for(m in 1:nrow(res_kegg)){
              str <- do.call(c,strsplit(res_kegg$geneID[m],split = "/"))
              str <- gene.df$SYMBOL[match(str,gene.df$ENTREZID)]
              str <- paste(str,collapse = "/")
              res_kegg$geneID[m] <- str
            }
            ench <- rbind(res_go,res_kegg)
          }else{
            ench <- res_go
          }
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
  data$cluster_new <- meta$cluster_new[match(colnames(data), substring(meta$file,4,21))]
  data$celltype <- "un_labeled"
  data$celltype <- meta$celltype[match(colnames(data), substring(meta$file,4,21))]
  data$sample <- sample2[i]
  
  dir.create(paste0(work_path, "/8_figure2/subcluster"), showWarnings = F)
  setwd(paste0(work_path, "/8_figure2/subcluster"))
  subdata <- data[,which(colnames(data) %in% meta$barcode)]
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
      if(is.null(kk)==F){
        res_kegg <- as.data.frame(kk@result)
        res_kegg <- res_kegg[which(res_kegg$pvalue<0.05),]
        for(m in 1:nrow(res_kegg)){
          str <- do.call(c,strsplit(res_kegg$geneID[m],split = "/"))
          str <- gene.df$SYMBOL[match(str,gene.df$ENTREZID)]
          str <- paste(str,collapse = "/")
          res_kegg$geneID[m] <- str
        }
        ench <- rbind(res_go,res_kegg)
      }else{
        ench <- res_go
      }
      ench$cluster <- cluster[k]
      rownames(ench) <- NULL
      list[[k]] <- ench
    }
  }
  ench_all <- do.call(rbind,list)
  write.csv(ench_all, "enrichment.csv", row.names = F)
}
















