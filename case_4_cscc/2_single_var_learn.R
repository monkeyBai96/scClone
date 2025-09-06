# 10000 vcf.gz
# ST
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
for(i in 1:2){
  ################################## path ########################################
  path = "H:/scclone/PMID32579974_ST"
  rna_path = paste0(path,"/paper/ST/GSE144239_RAW/",sample2[i])
  mut_path = paste0(path,"/",sample2[i],"/annovar")
  work_path = paste0(path,"/",sample2[i],"/work")
  setwd(paste0(work_path,"/2_svm"))
  ################################################################################
  
  
  ############################# SVM & signature ##################################
  # read data
  meta <- read.csv(paste0(work_path, "/0_exp/meta.csv"), row.names = 1)
  single_raw <- read.csv(paste0(work_path,"/1_basic/single_mut_vaf.csv"))
  #
  single_raw$var <- paste0(single_raw$Chr,"_",single_raw$Pos,"_",single_raw$Ref,">",single_raw$Alt)
  single_raw$DP <- ifelse(!is.na(single_raw$DP),single_raw$DP,single_raw$DPI)
  single_raw <- single_raw %>% filter(DP>10)
  single_snp <- single_raw[,c(10:36)]
  single_raw <- single_raw[,c("Chr","Pos","Ref","Alt","Qual","GQ","GQX","DP","DPI","AD","ADF","ADR","cell",
                              "Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","var")]
  single_raw$sample <- meta$sample[match(single_raw$cell, meta$file)]
  single_raw$cell <- meta$cell[match(single_raw$cell, meta$file)]
  
  # read RNA-editing
  pos_edit <- readRDS("H:/scclone/RNAedit/pos_edit_hg38.rds")
  single_snp$edit <- ifelse(single_raw$var %in% pos_edit$var, "yes", "no")
  extract_numeric <- function(x) {
    numeric_indices <- grepl("^\\s*[-+]?\\d*\\.?\\d+(e[-+]?\\d+)?\\s*$", x)
    if(all(numeric_indices==F)){ NULL }else{ x[numeric_indices] }
  }
  single_snp$dbsnp <- apply(single_snp[,1:27], 1, function(row) {
    numeric_elements <- extract_numeric(row)
    if(length(numeric_elements)>0){ max(numeric_elements) }else{ "." }
  })
  single_snp$dbsnp1 <- ifelse(as.numeric(single_snp$dbsnp)>0.5 & single_snp$dbsnp!=".","yes","no")
  
  # delete cell label
  single_svm <- single_raw[,c(18,3:12,16)]
  single_svm <- split_VAF_col(single_svm, col = "AD")
  single_svm <- split_VAF_col(single_svm, col = "ADF")
  single_svm <- split_VAF_col(single_svm, col = "ADR")
  single_svm$Ref <- ifelse(single_svm$Ref %in% c("A","T","C","G","-"), single_svm$Ref, "NN")
  single_svm$Alt <- ifelse(single_svm$Alt %in% c("A","T","C","G","-"), single_svm$Alt, "NN")
  single_svm$sbs <- paste0(single_svm$Ref,">",single_svm$Alt)
  length(unique(single_svm$var))
  single_svm <- single_svm[,c("var","sbs","Qual","GQ","GQX","DP","AD_1","AD_2","ADF_1","ADF_2","ADR_1","ADR_2","ExonicFunc.refGene")] # add ref & alt / sbs
  
  ################### positive:dbSNP;  negative:RNAedit DBS ######################
  table <- as.data.frame(table(single_raw$var))
  single_svm$mark <- ifelse(single_snp$edit=="no" & single_snp$dbsnp1=="yes" & single_raw$var %in% table$Var1[which(table$Freq>10)], "True", "Test")
  single_svm$mark <- ifelse(single_snp$edit=="yes" & single_snp$dbsnp1=="no", "False", single_svm$mark)
  single_svm$mark <- ifelse(nchar(single_raw$Ref)==nchar(single_raw$Alt) & nchar(single_raw$Alt)>1, "False", single_svm$mark)
  differences <- diff(single_raw$Pos)
  id_dbs <- sort(unique(c(which(differences==1),which(differences==1)+1)))
  id_snv <- which(single_raw$Ref %in% c("A","T","C","G") & single_raw$Alt %in% c("A","T","C","G"))
  id_dbs <- intersect(id_dbs, id_snv)
  single_svm$mark[id_dbs] <- "False"
  table(single_svm$mark)
  # prepare
  single_svm <- change_class(single_svm[,-1], factor_col = c("ExonicFunc.refGene","mark","sbs"))
  str(single_svm)
  rows_without_na <- apply(single_svm, 1, function(row) all(complete.cases(row)))
  id_complete <- which(rows_without_na==T)
  single_svm <- single_svm[id_complete,]
  single_snp <- single_snp[id_complete,]
  single_raw <- single_raw[id_complete,]
  
  # run SVM
  train <- single_svm[which(single_svm$mark %in% c("True","False")),]
  train$mark <- factor(train$mark, levels = c("True","False"))
  str(train)
  test <- single_svm[which(single_svm$mark %in% c("Test")),]
  str(test)
  # Fitting model
  fit <- svm(mark ~ ., data = train)
  # Predict Output
  pred <- predict(fit, test[,c(-13)])
  test$mark <- pred
  single_svm$mark1 <- ifelse(single_svm$mark=="Test", "Test", "Train")
  single_svm$mark1[which(single_svm$mark=="Test")] <- as.character(test$mark)
  table(single_svm$mark)
  table(single_svm$mark1)
  
  signature_T <- snv_to_sig(dat = single_raw[which(single_svm$mark1=="True" | single_svm$mark=="True"),c(1,2,3,4)], ref = "hg38") #
  signature_F <- snv_to_sig(dat = single_raw[which(single_svm$mark1=="False" | single_svm$mark=="False"),c(1,2,3,4)], ref = "hg38") #
  signature_A <- snv_to_sig(dat = single_raw[,c(1,2,3,4)], ref = "hg38") #
  
  # save
  p <- show_signature(array = as.array(table(signature_T$signature)), name = "SVM")
  ggsave("signature_T.pdf",p$plot, width = 8, height = 3)
  p <- show_signature(array = as.array(table(signature_F$signature)), name = "SVM")
  ggsave("signature_F.pdf",p$plot, width = 8, height = 3)
  p <- show_signature(array = as.array(table(signature_A$signature)), name = "SVM")
  ggsave("signature_A.pdf",p$plot, width = 8, height = 3)
  save <- cbind(single_raw[,c(1:4,13:18)], single_svm[,-which(colnames(single_svm)=="ExonicFunc.refGene")])
  write.csv(save, "single_svm.csv", row.names = F)
  write.csv(signature_T, "signature.csv", row.names = F)
}







