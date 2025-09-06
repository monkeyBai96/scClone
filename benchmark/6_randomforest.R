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
################################## path ########################################
path = "H:/scclone/PMID31558476_bm"
rna_path = paste0(path,"/paper")
mut_path = paste0(path,"/annovar")
work_path = paste0(path,"/work")
setwd(paste0(work_path,"/7_rf"))
################################################################################


# read data
meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
matrix_VAF <- read.csv(paste0(work_path, "/6_rpca/matrix_rpca.csv"), row.names = 1)
matrix_VAF <- rounding_VAF(matrix_VAF)

# run randomforest
demo <- run_randomforest(matrix = matrix_VAF, meta = meta, train_ratio = 0.7, output_path = paste0(work_path,"/7_rf"))
saveRDS(demo,"RF_result.rds")
importance <- demo[["importance"]]
var_keep <- rownames(importance)[which(importance$MeanDecreaseAccuracy>0&importance$MeanDecreaseGini>0)]
matrix_rf <- matrix_VAF[var_keep,]
meta_order <- meta %>% arrange(tissue, sample)
matrix_rf_order <- matrix_rf[,meta_order$cell]

# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:13)]
show_col(color)
names(color) <- c(unique(meta$sample),unique(meta$tissue))
anno_color <- list(sample = color[1:8], tissue = color[9:11])
graphics.off()
p1 <- pheatmap(matrix_rf, annotation_col = meta[,c("tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50),annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
p2 <- pheatmap(matrix_rf_order, annotation_col = meta_order[,c("tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50),annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = F)
p <- plot_grid(p1$gtable, p2$gtable, ncol = 2)
ggsave("heatmap_rf.pdf",p, width = 20, height = 10)
write.csv(matrix_rf, "matrix_rf.csv", row.names = T)

# plot distance
meta_order <- meta %>% arrange(tissue, sample)
matrix_order <- matrix_rf[,rownames(meta_order)]
distance <- as.matrix(suppressMessages(dist(t(matrix_order), method = "manhattan", diag = FALSE, upper = FALSE)))
order <- c()
group <- unique(meta_order$sample)
for(i in 1:length(group)){
  d <- dist(t(matrix_order[,which(meta_order$sample==group[i])]), method = "manhattan")
  x <- hclust(d, method = "complete")
  order <- c(order,x$labels[x$order])
}
p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("tissue","sample")], annotation_colors = anno_color,
                        annotation_row = meta_order[,c("tissue","sample")],
                        show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                        color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
ggsave("distance.pdf",p, width = 10, height = 10)








######################## rf training by sample #################################
dir.create(paste0(work_path,"/7_rf/by_sample"))
setwd(paste0(work_path,"/7_rf/by_sample"))

# read data
meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
matrix_VAF <- read.csv(paste0(work_path, "/6_rpca/matrix_rpca.csv"), row.names = 1)
matrix_VAF <- rounding_VAF(matrix_VAF)
meta$celltype <- meta$sample ###

# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:13)]
show_col(color)
names(color) <- c(unique(meta$sample),unique(meta$tissue))

# run randomforest
demo <- run_randomforest(matrix = matrix_VAF, meta = meta, train_ratio = 0.7, output_path = paste0(work_path,"/7_rf/by_sample"), celltype_color = color[1:8])
saveRDS(demo,"RF_result.rds")
importance <- demo[["importance"]]
var_keep <- rownames(importance)[which(importance$MeanDecreaseAccuracy>0.3 &importance$MeanDecreaseGini>0.3)]
matrix_rf <- matrix_VAF[var_keep,]
meta_order <- meta %>% arrange(tissue, sample)
matrix_rf_order <- matrix_rf[,meta_order$cell]

# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:13)]
show_col(color)
names(color) <- c(unique(meta$sample),unique(meta$tissue))
anno_color <- list(sample = color[1:8], tissue = color[9:11])
graphics.off()
p1 <- pheatmap(matrix_rf, annotation_col = meta[,c("tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50),annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
p2 <- pheatmap(matrix_rf_order, annotation_col = meta_order[,c("tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50),annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = F)
p <- plot_grid(p1$gtable, p2$gtable, ncol = 2)
ggsave("heatmap_rf.pdf",p, width = 20, height = 10)
write.csv(matrix_rf, "matrix_rf.csv", row.names = T)

# plot distance
meta_order <- meta %>% arrange(tissue, sample)
matrix_order <- matrix_rf[,rownames(meta_order)]
distance <- as.matrix(suppressMessages(dist(t(matrix_order), method = "manhattan", diag = FALSE, upper = FALSE)))
order <- c()
group <- unique(meta_order$sample)
for(i in 1:length(group)){
  d <- dist(t(matrix_order[,which(meta_order$sample==group[i])]), method = "manhattan")
  x <- hclust(d, method = "complete")
  order <- c(order,x$labels[x$order])
}
p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("tissue","sample")], annotation_colors = anno_color,
                        annotation_row = meta_order[,c("tissue","sample")],
                        show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                        color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
ggsave("distance.pdf",p, width = 10, height = 10)











