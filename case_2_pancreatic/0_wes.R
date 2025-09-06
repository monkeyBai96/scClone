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

# bulk
setwd("H:/scclone/PMID37612267/paper/wes")
data <- read.table("P3-PT_paired.hg19_multianno.txt", header = T ,sep = "\t", quote = "")
signature <- snv_to_sig(dat = data[,c(1,2,4,5)], ref = "hg19")
p1 <- show_signature(array = as.array(table(signature$signature)), name = "P3-PT-bulk")
setwd("H:/scclone/PMID37612267/paper/wes")
data <- read.table("P3-HM_paired.hg19_multianno.txt", header = T ,sep = "\t", quote = "")
signature <- snv_to_sig(dat = data[,c(1,2,4,5)], ref = "hg19")
p2 <- show_signature(array = as.array(table(signature$signature)), name = "P3-HM-bulk")
p <- cowplot::plot_grid(p1$plot,p2$plot,ncol = 1)
ggsave("bulk_signature.pdf",p, width = 9, height = 6)
signature <- cbind(p1$table, p2$table)
colnames(signature) <- c("P3-PT","P3-HM")
write.csv(signature, "bulk_signature.csv", row.names = T)

# single
single1 <- read.csv(paste0(path,"/",sample2[1],"/work/2_svm/signature.csv"))
p1 <- show_signature(array = as.array(table(single1$signature)), name = "SVM")
single2 <- read.csv(paste0(path,"/",sample2[2],"/work/2_svm/signature.csv"))
p2 <- show_signature(array = as.array(table(single2$signature)), name = "SVM")
single3 <- read.csv(paste0(path,"/",sample2[1],"/work/2_svm/signature_raw.csv"))
p3 <- show_signature(array = as.array(table(single3$signature)), name = "SVM")
single4 <- read.csv(paste0(path,"/",sample2[2],"/work/2_svm/signature_raw.csv"))
p4 <- show_signature(array = as.array(table(single4$signature)), name = "SVM")

signature <- cbind(p1$table, p2$table, p3$table, p4$table)
colnames(signature) <- c("P3_YF","P3_ZY","P3_YF_raw","P3_ZY_raw")
write.csv(signature, "single_signature.csv", row.names = T)


# compare
library(MutationalPatterns)
sig1 <- read.csv("bulk_signature.csv", row.names = 1)
sig2 <- read.csv("single_signature.csv", row.names = 1)
cos_sim <- cos_sim_matrix(sig1, sig2)
pdf("cos_sim.pdf",width = 4, height = 3)
plot_cosine_heatmap(cos_sim, cluster_rows = F, cluster_cols = F)
dev.off()
write.csv(cos_sim, "cos_sim.csv", row.names = T)











