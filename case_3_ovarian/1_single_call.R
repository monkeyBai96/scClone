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
sample7 <- c("P1","P2","P3","P4","P5","P6","P8")
for(i in 1:7){
  ################################################################################
  path = "H:/scclone/PMID38570491_ST"
  rna_path = paste0(path,"/paper/ST/",sample7[i])
  mut_path = paste0(path,"/",sample7[i],"/annovar")
  work_path = paste0(path,"/",sample7[i],"/work")
  setwd(paste0(work_path,"/1_basic"))
  ################################################################################
  
  
  ##################### read mutation data, per cell ###############################
  meta <- read.csv(paste0(work_path, "/0_exp/meta.csv"), row.names = 1)
  col_name <- c("Chr","Pos","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene",
                "ExonicFunc.refGene","AAChange.refGene","Qual","Filter","Info","Format","VAF","cell")
  data <- read_annovar(path = mut_path, col_keep = c(1,2,4,5,6,7,8,9,10,11:37,46:50))
  colnames(data)[c(1:9,37:42)] <- col_name
  data <- data %>% filter(cell %in% meta$file)
  sta <- data[,c(1:6,8,12,42)]
  data_VAF <- split_FORMAT_col(dat = data, split = 41, name = 40, sep = ":") # slow
  data_all <- cbind(data[,c(-40,-41)], data_VAF)
  
  # save
  write.csv(sta, "sta.csv", row.names = F)
  write.csv(data_all, "single_mut_vaf.csv", row.names = F)
  
  # statistic
  p <- statistic_single_mut(sta = sta, ref = "hg38")
  ggsave("basic_statistic.pdf", p, width = 23, height = 10)
}






