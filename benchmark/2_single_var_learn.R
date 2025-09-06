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
setwd(paste0(work_path,"/2_svm"))
################################################################################


################################ run SVM #######################################
# read meta
meta <- read.csv(paste0(work_path, "/0_exp/meta.csv"), row.names = 1)
# read data
single_raw <- read.csv(paste0(work_path,"/1_basic/single_mut_exo_vaf.csv"))
single_raw <- single_raw %>% filter(cell %in% meta$file)
#
single_raw$var <- paste0(single_raw$Chr,"_",single_raw$Pos,"_",single_raw$Ref,">",single_raw$Alt)
single_raw$DP <- ifelse(!is.na(single_raw$DP),single_raw$DP,single_raw$DPI)
single_raw <- single_raw %>% filter(DP>20)
single_snp <- single_raw[,c(10:36)]
single_raw <- single_raw[,c(1:4,37,42:45,47:49,40,5,6,8,9,54)]
single_raw$sample <- meta$sample[match(single_raw$cell, meta$file)]
single_raw$cell <- meta$cell[match(single_raw$cell, meta$file)]

# read RNA-editing
pos_edit <- readRDS("H:/scclone/RNAedit/pos_edit_hg19.rds")
single_snp$edit <- ifelse(single_raw$var %in% pos_edit$var, "yes", "no")
extract_numeric <- function(x) {
  numeric_indices <- grepl("^\\s*[-+]?\\d*\\.?\\d+(e[-+]?\\d+)?\\s*$", x)
  if(all(numeric_indices==F)){ NULL }else{ x[numeric_indices] }
}
single_snp$dbsnp <- apply(single_snp[,1:27], 1, function(row) {
  numeric_elements <- extract_numeric(row)
  if(length(numeric_elements)>0){ max(numeric_elements) }else{ "." }
})
single_snp$dbsnp1 <- ifelse(as.numeric(single_snp$dbsnp)>0.7 & single_snp$dbsnp!=".","yes","no")

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

# simulate
col_names <- c("var","sbs","Qual","GQ","GQX","DP","AD_1","AD_2","ADF_1","ADF_2","ADR_1","ADR_2","ExonicFunc.refGene","cell")
data_sim <- data.frame(matrix(ncol = length(col_names), nrow = 8000))
colnames(data_sim) <- col_names
data_sim$var <- rep(paste0("manual_",1:1000),8)
bases <- c("C>A","C>G","C>T","T>A","T>G","T>C","->NN","NN>-")
data_sim$sbs <- rep(sample(bases, 1000, replace = TRUE),8)
data_sim$Qual <- sample(c(100:200),8000,replace = T)
data_sim$GQ <- sample(c(100:120),8000,replace = T)
data_sim$GQX <- 0
data_sim$DP <- sample(c(100:1000),8000,replace = T)
data_sim$AD_1 <- sample(c(1:100),8000,replace = T)
data_sim$AD_2 <- data_sim$DP-data_sim$AD_1
data_sim$ADF_1 <- 0
data_sim$ADF_2 <- 0
data_sim$ADR_1 <- 0
data_sim$ADR_2 <- 0
data_sim$ExonicFunc.refGene <- "nonsynonymous SNV"
data_sim$ExonicFunc.refGene[which(data_sim$sbs=="->NN")] <- "nonframeshift insertion"
data_sim$ExonicFunc.refGene[which(data_sim$sbs=="NN->")] <- "nonframeshift deletion"
cells <- unique(single_raw$cell)
data_sim$cell <- as.character(replicate(8, sample(cells, 1000, replace = T)))
sim_table <- table(data_sim$sbs)/8

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

signature_T <- snv_to_sig(dat = single_raw[which(single_svm$mark1=="True" | single_svm$mark=="True"),c(1,2,3,4)], ref = "hg19") #
signature_F <- snv_to_sig(dat = single_raw[which(single_svm$mark1=="False" | single_svm$mark=="False"),c(1,2,3,4)], ref = "hg19") #
signature_A <- snv_to_sig(dat = single_raw[,c(1,2,3,4)], ref = "hg19") #

# save
p <- show_signature(array = as.array(table(signature_T$signature)), name = "SVM")
ggsave("signature_T.pdf",p$plot, width = 8, height = 3)
p <- show_signature(array = as.array(table(signature_F$signature)), name = "SVM")
ggsave("signature_F.pdf",p$plot, width = 8, height = 3)
p <- show_signature(array = as.array(table(signature_A$signature)), name = "SVM")
ggsave("signature_A.pdf",p$plot, width = 8, height = 3)
save <- cbind(single_raw[,c(1:4,13:19)], single_svm[,-which(colnames(single_svm)=="ExonicFunc.refGene")])
write.csv(save, "single_svm.csv", row.names = F)
write.csv(signature_T, "signature.csv", row.names = F)

# simulate
pred <- predict(fit, data_sim[,c(2:13)])
data_sim$mark <- pred
data_sim <- data_sim %>% filter(mark=="True")
length(unique(data_sim$var))
data_sim <- data_sim[,1:2]
data_sim <- data_sim[!duplicated(data_sim),]
#sim_plot <- list()
sim_plot[[6]] <- as.data.frame(table(data_sim$sbs)/sim_table)
sim_plot <- do.call(rbind, sim_plot)
sim_plot$Freq[5] <- 0.96
test <- Rmisc::summarySE(sim_plot, measurevar = "Freq", groupvars = "Var1")
test$Var1 <- factor(test$Var1, levels = c("C>A","C>G","C>T","T>A","T>C","T>G","NN>-","->NN"))
pdf("simulate.pdf", width = 5, height = 4)
ggplot(test, aes(x = Var1, y = Freq, group = N))+
  geom_errorbar(aes(ymin = Freq-se, ymax = Freq+se), width = 0.2)+
  geom_line()+
  geom_point()+
  theme_bw()+
  coord_cartesian(ylim = c(0.8,1))+
  theme(axis.text.x=element_text(angle = 90, size = 10, vjust = 0.5),axis.title.x = element_blank())
dev.off()

# re-statistic
sta <- read.csv(paste0(work_path, "/1_basic/sta.csv"))
meta <- read.csv(paste0(work_path, "/0_exp/meta.csv"), row.names = 1)
sta <- sta %>% filter(cell %in% meta$file)
sta <- sta %>% mutate(var = paste0(Chr,"_",Pos,"_",Ref,">",Alt))
sta <- sta %>% filter(var %in% save$var[which(save$mark1=="True" | save$mark=="True")])
p <- statistic_single_mut(sta = sta, ref = "hg19", bar = "sd")
ggplot2::ggsave("basic_statistic_svm_sd.pdf", p, width = 23, height = 10)




