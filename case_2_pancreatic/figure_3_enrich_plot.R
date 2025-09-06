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


i = 1
celltype = "Ductal_cell"
################################## path ########################################
rna_path = paste0(path,"/paper/single")
mut_path = paste0(path,"/",sample2[i],"/annovar")
work_path = paste0(path,"/",sample2[i],"/work")
################################################################################

setwd(paste0(work_path, "/8_figure"))
meta <- read.csv("meta_final.csv", row.names = 1)

# color
celltypes <- c("B_cell","Ductal_cell","Fibroblasts","Mac/Mono","Mast_cell","MKI67+cell","NK","Plasma_cell","T_cell")
celltypes_color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(1:9)]
names(celltypes_color) <- celltypes
celltype_color <- c(celltypes_color, "un_labeled"="gray")
cluster_new_color <- get_color_cluster_new(meta = meta, celltype_color = celltype_color)
cluster_new_color <- c(cluster_new_color, "un_labeled"="gray")

# plot
setwd(paste0(work_path, "/8_figure/subcluster/",celltype))
data <- read_excel("plot.xlsx", sheet = 1)
data <- data %>% arrange(cluster, Count)
data$Description <- factor(data$Description, levels = data$Description)

pdf("plot.pdf", width = 10, height = 5)
ggplot(data, aes(x = Description, y = log2(Count+1), fill = cluster))+
  geom_bar(stat = "identity", width = 0.8) + 
  xlab("") + ylab("log2 Gene count")+
  coord_flip() + 
  scale_fill_manual(values = cluster_new_color)+
  theme_bw(base_size = 15)
dev.off()
