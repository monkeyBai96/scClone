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
setwd(paste0(work_path, "/8_figure"))
################################################################################


index1 <- read.csv(paste0(work_path, "/4_vaf/index_final.csv"))
index2 <- read.csv(paste0(work_path, "/5_smooth/index_final.csv"))
index3 <- read.csv(paste0(work_path, "/6_rpca/index_final.csv"))
index <- rbind(melt(index1) %>% mutate(process = "vaf"),
               melt(index2) %>% mutate(process = "smooth"),
               melt(index3) %>% mutate(process = "rpca"))
index <- index[!is.na(index$value),]
colnames(index) <- c("sample","index","value","process")
index$process <- factor(index$process, levels = c("vaf","smooth","rpca"))
color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:15)]
show_col(color)
index$sample <- ifelse(index$sample=="all",as.character(index$index),index$sample)
names(color) <- unique(index$sample)

pdf("index1.pdf", width = 6, height = 4)
ggplot(index %>% filter(index=="si"), aes(x = sample, y = value, group = process, fill = sample, alpha = process))+
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8)+
  scale_fill_manual(values = color)+
  scale_alpha_manual(values = c(0.4,0.7,1))+
  facet_wrap(index ~ ., scales = "free", ncol = 4)+
  theme_bw(base_size = 12)
dev.off()
pdf("index2.pdf", width = 7, height = 4)
ggplot(index %>% filter(index!="si"), aes(x = sample, y = value, group = process, fill = sample, alpha = process))+
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8)+
  scale_fill_manual(values = color)+
  scale_alpha_manual(values = c(0.4,0.7,1))+
  facet_wrap(index ~ ., scales = "free", ncol = 5)+
  theme_bw(base_size = 12)
dev.off()




#cutoffs <- c(200,250,300,350,400,450,500,550,600)
#index1 <- list()
#for(i in 1:length(cutoffs)){
#  x <- read.csv(paste0(work_path, "/4_vaf/index_",cutoffs[i],".csv"))
#  x <- melt(x) %>% mutate(index = cutoffs[i], process = "vaf")
#  index1[[i]] <- x
#}
#index1 <- do.call(rbind, index1)
#index1 <- index1[!is.na(index1$value),]
#
#index2 <- list()
#for(i in 1:length(cutoffs)){
#  x <- read.csv(paste0(work_path, "/5_smooth/index_",cutoffs[i],".csv"))
#  x <- melt(x) %>% mutate(index = cutoffs[i], process = "smooth")
#  index2[[i]] <- x
#}
#index2 <- do.call(rbind, index2)
#index2 <- index2[!is.na(index2$value),]
#
#index3 <- list()
#for(i in 1:length(cutoffs)){
#  x <- read.csv(paste0(work_path, "/6_rpca/index_",cutoffs[i],".csv"))
#  x <- melt(x) %>% mutate(index = cutoffs[i], process = "rpca")
#  index3[[i]] <- x
#}
#index3 <- do.call(rbind, index3)
#index3 <- index3[!is.na(index3$value),]
#index <- rbind(index1, index2, index3)
#
#color <- c(pal_npg("nrc")(10), pal_nejm("default")(10))[c(5,7,6,3,1,4,2,9, 11:14)]
#show_col(color)
#index$sample <- ifelse(index$sample=="all",as.character(index$variable),index$sample)
#names(color) <- unique(index$sample)
#index$process <- factor(index$process, levels = c("vaf","smooth","rpca"))
#ggplot(index %>% filter(variable=="si"), aes(x = process, y = value, group = index, fill = sample, alpha = process))+
#  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8)+
#  scale_fill_manual(values = color)+
#  scale_alpha_manual(values = c(0.4,0.7,1))+
#  facet_wrap(sample ~ ., scales = "free", ncol = 4)+
#  theme_bw(base_size = 12)
#pdf("demo2025index.pdf", height = 20, width = 7)
#ggplot(index %>% filter(variable!="si"), aes(x = process, y = value, group = index, fill = sample, alpha = process))+
#  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8)+
#  scale_fill_manual(values = color)+
#  scale_alpha_manual(values = c(0.4,0.7,1))+
#  facet_wrap(index ~ variable, scales = "free", ncol = 4)+
#  theme_bw(base_size = 12)
#dev.off()
