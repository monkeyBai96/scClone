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


i=1
celltype = "type_1"
################################################################################
path = "H:/scclone/PMID32579974_ST"
rna_path = paste0(path,"/paper/ST/GSE144239_RAW/",sample2[i])
mut_path = paste0(path,"/",sample2[i],"/annovar")
work_path = paste0(path,"/",sample2[i],"/work")
################################################################################

setwd(paste0(work_path, "/8_figure"))
meta <- read.csv("meta_final.csv", row.names = 1)

# color
celltypes <- paste0("type_",c(0:11))
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(8))
show_col(celltypes_color)
celltypes_color <- celltypes_color[c(1:5,13,12,9,11,6,7,10)]
names(celltypes_color) <- celltypes
celltype_color <- c(celltypes_color, "un_labeled"="gray")
cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")

# plot
setwd(paste0(work_path, "/8_figure/subcluster/",celltype))
data <- read_excel("plot.xlsx", sheet = 1)
data <- data %>% arrange(Count, cluster)
data <- rbind(data[which(data$cluster=="type_1_2"),], data[which(data$cluster=="type_1_1"),])
data$Description <- factor(data$Description, levels = data$Description)

pdf("plot.pdf", width = 8, height = 5)
ggplot(data, aes(x = Description, y = log2(Count+1), fill = cluster))+
  geom_bar(stat = "identity", width = 0.8) + 
  xlab("") + ylab("log2 Gene count")+
  coord_flip() + 
  scale_fill_manual(values = cluster_new_color)+
  theme_bw(base_size = 15)
dev.off()














i=7
################################################################################
path = "H:/scclone/PMID38570491_ST"
rna_path = paste0(path,"/paper/ST/",sample7[i])
mut_path = paste0(path,"/",sample7[i],"/annovar")
work_path = paste0(path,"/",sample7[i],"/work")
################################################################################

setwd(paste0(work_path, "/8_figure2"))
meta <- read.csv("meta_final.csv", row.names = 1)

# color
celltypes <- paste0("type_",c(0:11))
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(8))
show_col(celltypes_color)
celltypes_color <- celltypes_color[c(1:5,13,12,9,11,6,7,10)]
names(celltypes_color) <- celltypes
celltype_color <- c(celltypes_color, "un_labeled"="gray")
cluster_new_color <- c("cluster_1"=pal_nejm("default")(8)[5], "cluster_2"=pal_nejm("default")(8)[7])
cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")

# plot
setwd(paste0(work_path, "/8_figure2/subcluster"))
data <- read_excel("plot.xlsx", sheet = 1)
data <- data %>% arrange(Count)
data <- rbind(data[which(data$cluster=="cluster_2"),], data[which(data$cluster=="cluster_1"),])
data$Description <- factor(data$Description, levels = data$Description)

pdf("plot.pdf", width = 8, height = 5)
ggplot(data, aes(x = Description, y = log2(Count+1), fill = cluster))+
  geom_bar(stat = "identity", width = 0.8) + 
  xlab("") + ylab("log2 Gene count")+
  coord_flip() + 
  scale_fill_manual(values = cluster_new_color)+
  theme_bw(base_size = 15)
dev.off()




