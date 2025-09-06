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
setwd(paste0(work_path,"/5_smooth"))
################################################################################

# read
meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)
matrix_VAF <- read.csv(paste0(work_path,"/4_vaf/matrix_VAF.csv"), row.names = 1)
matrix_ALT <- read.csv(paste0(work_path,"/4_vaf/matrix_ALT.csv"), row.names = 1)
matrix_REF <- read.csv(paste0(work_path,"/4_vaf/matrix_REF.csv"), row.names = 1)
matrix_VAF[is.na(matrix_VAF)] <- 0
matrix_VAF[matrix_VAF>2] <- 2

# hclust
demo <- as.matrix(dist(t(matrix_VAF), method = "manhattan", diag = FALSE, upper = FALSE))
hist(unlist(demo), xlim = c(0,max(demo)), breaks = 100)
cluster <- get_hclust(matrix = matrix_VAF, cut = 300, method = "manhattan")
meta$cluster_tem <- cluster[match(meta$cell, names(cluster))]
table(meta$cluster_tem)

# refine meta: celltype and cluster_tem
meta_order <- meta %>% arrange(tissue, sample)
meta_select <- refine_meta(meta = meta_order, min_cell = 5)
table(meta_select$cluster_tem, meta_select$sample)

# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9,14:15, 11:13)]
show_col(color)
names(color) <- c(unique(meta$sample),unique(meta$tissue))
anno_color <- list(sample = color[1:10], tissue = color[11:13])
p1 <- pheatmap(matrix_VAF, annotation_col = meta_order[,c("cluster_tem","tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
p2 <- pheatmap(matrix_VAF[,meta_order$cell], annotation_col = meta_order[,c("cluster_tem","tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = F)

# signature
signature <- read.csv(paste0(work_path,"/2_svm/signature.csv"))
sig <- table(signature$signature)

# slow
matrix_VAF <- matrix_VAF[,which(colnames(matrix_VAF) %in% rownames(meta_select))]
matrix_ALT <- matrix_ALT[,which(colnames(matrix_ALT) %in% rownames(meta_select))]
matrix_REF <- matrix_REF[,which(colnames(matrix_REF) %in% rownames(meta_select))]
matrix_smooth <- make_smooth_matrix(matrix_VAF = matrix_VAF, matrix_ALT = matrix_ALT,matrix_REF = matrix_REF,meta = meta_select, sig = sig,
                                    ref = "hg19", 
                                    method = "manhattan",
                                    leader_min_depth = 20,
                                    neighbor_min_dis = 300,
                                    K_smooth = 0.2)
graphics.off()
p3 <- pheatmap(matrix_smooth, annotation_col = meta_select[,c("cluster_tem","tissue","sample")], color = colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
                show_rownames = F, show_colnames = F, cluster_cols = T, clustering_distance_cols = "manhattan")
p4 <- pheatmap(matrix_smooth[,meta_select$cell], annotation_col = meta_select[,c("cluster_tem","tissue","sample")], colorRampPalette(brewer.pal(9,"GnBu"))(50), annotation_colors = anno_color,
               show_rownames = F, show_colnames = F, cluster_cols = F)

# save
write.csv(matrix_smooth, "matrix_smooth.csv", row.names = T)
write.csv(matrix_VAF, "matrix_VAF.csv", row.names = T)
write.csv(matrix_ALT, "matrix_ALT.csv", row.names = T)
write.csv(matrix_REF, "matrix_REF.csv", row.names = T)
write.csv(meta_select, paste0(work_path,"/1_basic/meta_select.csv"), row.names = T)
p <- plot_grid(p1$gtable, p2$gtable, p3$gtable, p4$gtable, ncol = 2)
ggsave("heatmap_smooth.pdf", p, width = 20, height = 20)






# plot VAF calculate
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
setwd(paste0(work_path,"/5_smooth"))
################################################################################


matrix_ALT <- read.csv("matrix_ALT.csv", row.names = 1)
matrix_REF <- read.csv("matrix_REF.csv", row.names = 1)
matrix_smooth <- read.csv("matrix_smooth.csv", row.names = 1)
single <- read.csv(paste0(work_path, "/3_filter/single_filter_more.csv"))
meta <- read.csv(paste0(work_path,"/1_basic/meta_select.csv"), row.names = 1)

single$VAF <- single$AD_2/(single$AD_1+single$AD_2)
matrix_VAF_raw <- col1_to_rowname(as.matrix(dcast(single[,c("cell","var","VAF")], var~cell, value.var = "VAF", fun.aggregate=sum)))
matrix_VAF_raw[is.na(matrix_VAF_raw)] <- 0
matrix_VAF_raw <- matrix_VAF_raw[rownames(matrix_smooth),colnames(matrix_smooth)]
data_tem <- data.frame(VAF = as.numeric(unlist(matrix_VAF_raw)),
                       ALT = as.numeric(unlist(matrix_ALT)),
                       REF = as.numeric(unlist(matrix_REF)),
                       genotype = as.numeric(unlist(matrix_smooth)))
data_tem <- data_tem %>% mutate(alt_depth = ifelse(ALT==0, "drop-out",ifelse(ALT+REF<10,"<10",">10")))
data_tem$alt_depth <- factor(data_tem$alt_depth, levels = c(">10","<10","drop-out"))
color <- pal_npg("nrc")(10)[c(2,1,3)]
p1 <- ggplot(data_tem, aes(VAF, fill = alt_depth, group = alt_depth))+
  geom_histogram(position = "stack")+
  scale_fill_manual(values = color)+
  theme_few()
p2 <- ggplot(data_tem, aes(genotype, fill = alt_depth, group = alt_depth))+
  geom_histogram(position = "stack")+
  scale_fill_manual(values = color)+
  theme_few()
graphics.off()
pdf("hist_VAF_raw.pdf", width = 6, height = 5)
gg.gap(plot=p1,segments=list(c(10000,20000)),ylim=c(0,3*10^5))
#add.legend(plot = p1,margin = c(top=1,right=1,bottom=200,left=300))
dev.off()
pdf("hist_VAF.pdf", width = 6, height = 5)
gg.gap(plot=p2,segments=list(c(35000,40000)),ylim=c(0,3*10^5))
#add.legend(plot = p2,margin = c(top=1,right=1,bottom=200,left=300))
dev.off()
pdf("hist_VAF_new.pdf", width = 6, height = 5)
p2
dev.off()
# plot distance
meta_order <- meta %>% arrange(tissue, sample)
rownames(meta_order) <- meta_order$cell
matrix_order <- matrix_smooth[,rownames(meta_order)]
distance <- as.matrix(suppressMessages(dist(t(matrix_order), method = "manhattan", diag = FALSE, upper = FALSE)))
order <- c()
group <- unique(meta_order$sample)
for(j in 1:length(group)){
    d <- dist(t(matrix_order[,which(meta_order$sample==group[j])]), method = "manhattan", diag = FALSE, upper = FALSE)
    if(nrow(d)>2){
      x <- hclust(d, method = "complete")
      order <- c(order,x$labels[x$order])
    }else{
      order <- c(order,meta_order$cell[which(meta_order$sample==group[j])])
    }
  }
# color
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9,14:15, 11:13)]
show_col(color)
names(color) <- c(unique(meta$sample),unique(meta$tissue))
anno_color <- list(sample = color[1:10], tissue = color[11:13])
p <- pheatmap::pheatmap(distance[order,order], annotation_col = meta_order[,c("tissue","sample")], annotation_colors = anno_color,
                        annotation_row = meta_order[,c("tissue","sample")],
                        show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F,
                        color = colorRampPalette(brewer.pal(11,"RdBu"))(50))
ggsave("distance.pdf", p, width = 10, height = 10)

# statistic
demo <- as.matrix(dist(t(matrix_smooth), method = "manhattan", diag = FALSE, upper = FALSE))
hist(unlist(demo), xlim = c(0,max(demo)), breaks = 100)
cluster <- get_hclust(matrix = matrix_smooth, cut = 300, method = "manhattan")
table(cluster)
meta_order$cluster_tem <- cluster[match(meta_order$cell, names(cluster))]
index <- get_assess_index(meta_order = meta_order, matrix_order = matrix_order)
write.csv(index, paste0("index_final.csv"), row.names = F)
for(i in c(200,250,300,350,400,450,500,550,600)){
  cluster <- get_hclust(matrix = matrix_smooth, cut = i, method = "manhattan")
  table(cluster)
  
  meta_order$cluster_tem <- cluster[match(meta_order$cell, names(cluster))]
  index <- get_assess_index(meta_order = meta_order, matrix_order = matrix_order)
  write.csv(index, paste0("index_",i,".csv"), row.names = F)
}







