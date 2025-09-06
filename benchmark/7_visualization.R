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
setwd(paste0(work_path,"/8_figure"))
################################################################################


# read data
meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
matrix_VAF <- read.csv(paste0(work_path,"/7_rf/matrix_rf.csv"), row.names = 1)
matrix_VAF[is.na(matrix_VAF)] <- 0
matrix_VAF[matrix_VAF>2] <- 2

cluster <- get_hclust(matrix = matrix_VAF, method = "manhattan", cut = 300)
meta$cluster_new <- paste0("cluster_",cluster[match(meta$cell, names(cluster))])
table(meta$cluster_new ,meta$sample)
table <- as.data.frame(table(paste0(meta$cluster_new, "_",meta$sample)))
meta_select <- meta %>% filter(paste0(meta$cluster_new, "_",meta$sample) %in% table$Var1[which(table$Freq>2)]) %>% arrange(tissue, sample, cluster_new)


# color
#celltype_color <- get_color_celltype(meta = meta_select, pattern = "Set1")
#cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
#heatmap_color <- list(celltype = celltype_color, cluster_new = cluster_new_color)
#sankey_color <- c(celltype_color,cluster_new_color)

# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:13)]
names(color) <- c(unique(meta$sample),unique(meta$tissue))
show_col(color)
cluster_new_color <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta_select$cluster_new)))[c(1:length(unique(meta_select$cluster_new)))]
names(cluster_new_color) <- c(unique(meta$cluster_new))
show_col(cluster_new_color)
heatmap_color <- list(sample = color[1:8], tissue = color[9:11], cluster_new = cluster_new_color)
sankey_color <- c(color, cluster_new_color)

p1 <- pheatmap(matrix_VAF[,meta_select$cell], annotation_col = meta_select[,c("cluster_new","sample","tissue")], show_rownames = F, show_colnames = F, 
               cluster_cols = F, cluster_rows = F, color = colorRampPalette(brewer.pal(9,"GnBu"))(50),
               clustering_distance_cols = "manhattan", annotation_colors = heatmap_color)

df <- meta_select %>% make_long(tissue, sample, cluster_new) 
df$node <- factor(df$node, levels = c(unique(meta_select$tissue),
                                      unique(meta_select$sample),
                                      unique(meta_select$cluster_new)))
p2 <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.7, node.color = "gray70", smooth = 5, width = 0.2) +
  geom_sankey_label(aes(fill = factor(node)), size = 5, color = "white") +
  scale_fill_manual(values = sankey_color)+
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
  ggtitle("tissue-sample-cluster")
p2

table <- as.data.frame(table(meta_select$cluster_new, meta_select$sample))
p3 <- ggplot(table, aes(area = Freq, fill = Var1, label = Var1, subgroup = Var2)) +
  geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",grow = TRUE,alpha=.6)+
  geom_treemap_subgroup_border(color='white')+
  scale_fill_manual(values = sankey_color)
p3

# save
graphics.off()
p <- plot_grid(p1$gtable, p2, p3, ncol = 3)
ggsave("demo.pdf",p, width = 30, height = 10)
write.csv(meta_select, "meta_final.csv", row.names = T)
write.csv(matrix_VAF[,meta_select$cell], "matrix_final.csv")






######################## rf training by sample #################################
dir.create(paste0(work_path,"/8_figure/by_sample"))
setwd(paste0(work_path,"/8_figure/by_sample"))

# read data
meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
matrix_VAF <- read.csv(paste0(work_path,"/7_rf/by_sample/matrix_rf.csv"), row.names = 1)
matrix_VAF[is.na(matrix_VAF)] <- 0
matrix_VAF[matrix_VAF>2] <- 2

cluster <- get_hclust(matrix = matrix_VAF, method = "manhattan", cut = 90)
meta$cluster_new <- paste0("cluster_",cluster[match(meta$cell, names(cluster))])
table(meta$cluster_new ,meta$sample)
table <- as.data.frame(table(paste0(meta$cluster_new, "_",meta$sample)))
meta_select <- meta %>% filter(paste0(meta$cluster_new, "_",meta$sample) %in% table$Var1[which(table$Freq>3)]) %>% arrange(tissue, sample, cluster_new)

# color
#celltype_color <- get_color_celltype(meta = meta_select, pattern = "Set1")
#cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
#heatmap_color <- list(celltype = celltype_color, cluster_new = cluster_new_color)
#sankey_color <- c(celltype_color,cluster_new_color)

# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:13)]
names(color) <- c(unique(meta$sample),unique(meta$tissue))
show_col(color)
cluster_new_color <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta_select$cluster_new)))[c(1:length(unique(meta_select$cluster_new)))]
names(cluster_new_color) <- c(unique(meta_select$cluster_new))
show_col(cluster_new_color)
heatmap_color <- list(sample = color[1:8], tissue = color[9:11], cluster_new = cluster_new_color)
sankey_color <- c(color, cluster_new_color)

p1 <- pheatmap(matrix_VAF[,meta_select$cell], annotation_col = meta_select[,c("cluster_new","sample","tissue")], show_rownames = F, show_colnames = F, 
               cluster_cols = F, cluster_rows = T, color = colorRampPalette(brewer.pal(9,"GnBu"))(50),
               clustering_distance_cols = "manhattan", annotation_colors = heatmap_color)

p11 <- pheatmap(matrix_VAF[,meta_select$cell], annotation_col = meta_select[,c("cluster_new","sample","tissue")], show_rownames = F, show_colnames = F, 
               cluster_cols = T, cluster_rows = T, color = colorRampPalette(brewer.pal(9,"GnBu"))(50),
               clustering_distance_cols = "manhattan", annotation_colors = heatmap_color)

df <- meta_select %>% make_long(tissue, sample, cluster_new) 
df$node <- factor(df$node, levels = c(unique(meta_select$tissue),
                                      unique(meta_select$sample),
                                      unique(meta_select$cluster_new)))
p2 <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.7, node.color = "gray70", smooth = 5, width = 0.2) +
  geom_sankey_label(aes(fill = factor(node)), size = 5, color = "white") +
  scale_fill_manual(values = sankey_color)+
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
  ggtitle("tissue-sample-cluster")
p2

table <- as.data.frame(table(meta_select$cluster_new, meta_select$sample))
p3 <- ggplot(table, aes(area = Freq, fill = Var1, label = Var1, subgroup = Var2)) +
  geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",grow = TRUE,alpha=.6)+
  geom_treemap_subgroup_border(color='white')+
  scale_fill_manual(values = sankey_color)
p3

# save
graphics.off()
p <- plot_grid(p1$gtable,p11$gtable, p2, p3, ncol = 4)
ggsave("demo.pdf",p, width = 40, height = 10)
write.csv(meta_select, "meta_final.csv", row.names = T)
write.csv(matrix_VAF[,meta_select$cell], "matrix_final.csv")

mut <- data.frame(var = rownames(matrix_VAF))
mut <- separate(mut, var, into = c("chr","pos","refalt"), sep = "_")
mut <- separate(mut, refalt, into = c("ref","alt"), sep = ">")
signature <- snv_to_sig(dat = mut, ref = "hg19")
p <- show_signature(array = as.array(table(signature$signature)), name = "RF")
ggsave("signature_final.pdf",p$plot, width = 8, height = 3)







