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




i = 2
celltype = "Tumor_cell"
################################## path ########################################
rna_path = paste0(path,"/paper")
mut_path = paste0(path,"/",HCC_name[i],"/annovar")
work_path = paste0(path,"/",HCC_name[i],"/work")
################################################################################

setwd(paste0(work_path, "/8_figure"))
meta <- read.csv("meta_final.csv", row.names = 1)

# color
n = length(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell", "un_labeled"))
color <- c(pal_npg("nrc")(10))[c(1:n)]
show_col(color)
names(color) <- c(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell", "un_labeled"))
celltype_color <- color
cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")

# plot
setwd(paste0(work_path, "/8_figure/subcluster/",celltype))
data <- read_excel("plot.xlsx", sheet = 1)
data <- data %>% arrange(cluster, Count)
data$Description <- factor(data$Description, levels = data$Description)

pdf("plot.pdf", width = 8, height = 6)
ggplot(data, aes(x = Description, y = log2(Count+1), fill = cluster))+
  geom_bar(stat = "identity", width = 0.8) + 
  xlab("") + ylab("log2 Gene count")+
  coord_flip() + 
  scale_fill_manual(values = cluster_new_color)+
  theme_bw(base_size = 15)
dev.off()







i = 3
celltype = "T_cell"
################################## path ########################################
rna_path = paste0(path,"/paper")
mut_path = paste0(path,"/",HCC_name[i],"/annovar")
work_path = paste0(path,"/",HCC_name[i],"/work")
################################################################################

setwd(paste0(work_path, "/8_figure"))
meta <- read.csv("meta_final.csv", row.names = 1)

# color
n = length(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell", "un_labeled"))
color <- c(pal_npg("nrc")(10))[c(1:n)]
show_col(color)
names(color) <- c(c("Tumor_cell","Others_2","Macrophage","Others_10","T_cell","B_cell", "un_labeled"))
celltype_color <- color
cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")

# plot
setwd(paste0(work_path, "/8_figure/subcluster/",celltype))
data <- read_excel("plot.xlsx", sheet = 1)
data <- data %>% arrange(Count)
data$Description <- factor(data$Description, levels = data$Description)

pdf("plot.pdf", width = 6.5, height = 4)
ggplot(data, aes(x = Description, y = log2(Count+1), fill = cluster))+
  geom_bar(stat = "identity", width = 0.8) + 
  xlab("") + ylab("log2 Gene count")+
  coord_flip() + 
  scale_fill_manual(values = cluster_new_color)+
  theme_bw(base_size = 15)
dev.off()
