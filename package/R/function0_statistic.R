# functions for scClone
# bsh2024


#' @title Single-cell mutation statistics
#' @description Plot mutation types, substitutions, RNA editing, dbSNP and cumulative variants.
#' @param sta Mutation data frame.
#' @param ref Genome build hg19/hg38/mm10.
#' @param bar Error bar type se or sd.
#' @return ggplot object
#' @export
statistic_single_mut <- function(sta = sta, ref = ref, bar = "se"){
  if(bar=="se"){
    ebtop<-function(x){return(mean(x)+sd(x)/sqrt(length(x)))}
    ebbottom<-function(x){return(mean(x)-sd(x)/sqrt(length(x)))}
  }else{
    ebtop<-function(x){return(mean(x)+sd(x))}
    ebbottom<-function(x){return(mean(x)-sd(x))}
  }

  cat("Printing Func.refGene and substitutes...\n")
  # statistic Func.refGene
  table <- as.data.frame(table(sta$Func.refGene, sta$cell)) %>% mutate(Var1=as.character(Var1))
  table$Var1[which(grepl("ncRNA",table$Var1)==T)] <- "ncRNA"
  table$Var1[which(grepl(";",table$Var1)==T)] <- "others"
  table_sum <- table %>%
    group_by(Var1, Var2)  %>%
    summarise(Freq = sum(Freq, na.rm = TRUE))
  Func.refGene_order <- c("exonic","splicing","intronic","UTR3","UTR5","downstream","upstream","intergenic","ncRNA","others")
  table_sum$Var1 <- factor(table_sum$Var1, levels = Func.refGene_order)
  p1 <- ggplot(data = table_sum, aes(x=Var1, y=Freq, fill=Var1))+
    stat_summary(geom = "bar",fun = "mean",position = position_dodge(0.9))+
    stat_summary(geom = "errorbar",fun.min = ebbottom,fun.max = ebtop,position = position_dodge(0.9),width=0.5, size = 0.5)+
    scale_fill_npg()+
    theme_bw()+
    theme(axis.text.x=element_text(angle = 90, size = 15, vjust = 0.5),axis.title.x = element_blank())

  # statistic ATCG
  sta$atcg <- ifelse(sta$Ref%in%c("A","T","C","G")&sta$Alt%in%c("A","T","C","G"),paste0(sta$Ref,">",sta$Alt),
                     ifelse(sta$Ref=="-","INS", ifelse(sta$Alt=="-","DEL", "Others")))
  table <- as.data.frame(table(sta$atcg, sta$cell)) %>% mutate(Var1=as.character(Var1))
  table$Var1[which(table$Var1=="G>T")] <- "C>A"
  table$Var1[which(table$Var1=="G>C")] <- "C>G"
  table$Var1[which(table$Var1=="G>A")] <- "C>T"
  table$Var1[which(table$Var1=="A>T")] <- "T>A"
  table$Var1[which(table$Var1=="A>G")] <- "T>C"
  table$Var1[which(table$Var1=="A>C")] <- "T>G"
  table <- table[which(table$Var1%in%c("C>A","C>G","C>T","DEL","INS","T>A","T>C","T>G")),]
  table_sum <- table %>%
    group_by(Var1, Var2) %>%
    summarise(Freq = sum(Freq, na.rm = TRUE))
  ATCG_order <- c("C>A","C>G","C>T","T>A","T>C","T>G","INS","DEL")
  table_sum$Var1 <- factor(table_sum$Var1, levels = ATCG_order)
  p2 <- ggplot(data = table_sum, aes(x=Var1, y=Freq, fill=Var1))+
    stat_summary(geom = "bar",fun = "mean",position = position_dodge(0.9))+
    stat_summary(geom = "errorbar",fun.min = ebbottom,fun.max = ebtop,position = position_dodge(0.9),width=0.5, size = 0.5)+
    scale_fill_npg()+
    theme_bw()+
    theme(axis.text.x=element_text(angle = 90, size = 15, vjust = 0.5),axis.title.x = element_blank())

  # statistic RNAediting & dbSNP
  cat("Printing RNAediting dbSNP and depth...\n")
  if(ref=="hg19"){
    pos_edit <- readRDS("H:/scclone/RNAedit/pos_edit_hg19.rds")
  }else if(ref=="hg38"){
    pos_edit <- readRDS("H:/scclone/RNAedit/pos_edit_hg38.rds")
  }else if(ref=="mm10"){
    pos_edit <- readRDS("H:/scclone/RNAedit/pos_edit_mm10.rds")
  }else{
    cat("Error ref provided.\n")
  }
  sta$var <- paste0(sta$Chr,"_",sta$Pos,"_",sta$Ref,">",sta$Alt)
  sta$edit <- ifelse(sta$var %in% pos_edit$var, "Edit", "unknown")
  sta$dbsnp <- ifelse(sta$avsnp147!=".", "dbSNP", "unknown")
  table_var <- as.data.frame(table(sta$var)) %>% mutate(Var1=as.character(Var1))
  sta$label <- ifelse(sta$var %in% table_var$Var1[which(table_var$Freq>1)], ">1", 1)
  sta$label <- ifelse(sta$var %in% table_var$Var1[which(table_var$Freq>5)], ">5", sta$label)
  sta$label <- ifelse(sta$var %in% table_var$Var1[which(table_var$Freq>20)], ">20", sta$label)
  pie <- rbind(as.data.frame(table(sta$edit)), as.data.frame(table(sta$dbsnp)), as.data.frame(table(sta$label)))
  pie <- pie %>% mutate(prob = Freq/nrow(sta)*100) %>% mutate(group = c("edit","edit","dbsnp","dbsnp","label","label","label","label"))

  pie_tem <- pie %>% filter(group=="edit") %>%
    mutate(Var1 = factor(Var1, levels = c("Edit","unknown"))) %>%
    mutate(pos = 100-cumsum(prob)+0.5*prob)
  p3 <- ggplot(pie_tem, aes(x = "", y = prob, fill = Var1 )) +
    geom_bar(width = 1, stat = "identity", color = "white")+
    coord_polar("y", start = 0)+
    geom_text(aes(y = pos, label = signif(prob,3)), color = "black", size = 5)+
    scale_fill_npg() +
    theme_void()

  pie_tem <- pie %>% filter(group=="dbsnp") %>%
    mutate(Var1 = factor(Var1, levels = c("dbSNP","unknown"))) %>%
    mutate(pos = 100-cumsum(prob)+0.5*prob)
  p4 <- ggplot(pie_tem, aes(x = "", y = prob, fill = Var1 )) +
    geom_bar(width = 1, stat = "identity", color = "white")+
    coord_polar("y", start = 0)+
    geom_text(aes(y = pos, label = signif(prob,3)), color = "black", size = 5)+
    scale_fill_npg() +
    theme_void()

  pie_tem <- pie %>% filter(group=="label") %>%
    mutate(Var1 = factor(Var1, levels = c(">20",">5",">1","1"))) %>%
    mutate(pos = 100-cumsum(prob)+0.5*prob)
  p5 <- ggplot(pie_tem, aes(x = "", y = prob, fill = Var1 )) +
    geom_bar(width = 1, stat = "identity", color = "white")+
    coord_polar("y", start = 0)+
    geom_text(aes(y = pos, label = signif(prob,3)), color = "black", size = 5)+
    scale_fill_npg() +
    theme_void()

  # cumsum variants
  table_var <- sta %>%
    group_by(cell) %>%
    summarise(var = n_distinct(var)) %>%
    mutate(var_cumsum = NA, no = 1:length(unique(cell)))
  idx <- unique(c(seq(1,nrow(table_var),20), nrow(table_var)))
  pb <- txtProgressBar(min = 0, max = nrow(table_var), style = 3)
  var_cumsum <- c()
  for(i in 1:length(idx)){
    setTxtProgressBar(pb, idx[i])
    if(i!=1){
      tem <- sta %>% filter(cell %in% table_var$cell[(idx[i-1]+1):idx[i]])
    }else{
      tem <- sta %>% filter(cell %in% table_var$cell[1:idx[i]])
    }
    tem_var <- unique(tem$var)
    var_cumsum <- unique(c(var_cumsum,tem_var))
    table_var$var_cumsum[idx[i]] <- length(var_cumsum)
  }
  close(pb)
  scale <- ceiling(table_var$var_cumsum[nrow(table_var)]/max(table_var$var))

  p6 <- ggplot(table_var, aes(x = "", y = var ))+
    geom_boxplot(width=0.2, outlier.shape = NA, fill = pal_npg("nrc", alpha = 0.6)(4)[2])+
    geom_jitter(width=0.2, color=pal_npg("nrc", alpha=0.6)(4)[3], alpha = 0.8)+
    geom_line(data = table_var[idx,], aes(x=seq(from=1-0.5, to=1+0.5, length.out=nrow(table_var))[idx], y=var_cumsum/scale),
              color=pal_npg("nrc", alpha = 0.6)(4)[4], linetype="solid", linewidth = 3, alpha = 0.5)+
    theme_bw()+
    theme(axis.text.x=element_text(size = 15, vjust = 0.5),axis.title.x = element_blank(), axis.ticks.x = element_blank())+
    scale_y_continuous(expand = c(0.06,0), sec.axis = sec_axis(~.*scale, name = "cumsum"))

  p <- patchwork::wrap_plots(p1, p2, p6, p3, p4, p5, byrow = T, ncol = 3)
  return(p)
  cat("DONE! \n")
}


#' @title Cluster quality indices
#' @description Compute silhouette, CH, DB, ARI, FM and NMI for sample vs cluster labels.
#' @param meta_order Metadata with sample and cluster_tem.
#' @param matrix_order Expression matrix rows ordered by meta_order.
#' @return Data frame of indices per sample and overall.
#' @export
#' @importFrom fpc calinhara
#' @importFrom clusterSim index.DB
#' @importFrom mclust adjustedRandIndex
#' @importFrom dendextend FM_index
#' @importFrom aricode NMI
get_assess_index <- function(meta_order = meta_order, matrix_order = matrix_order){

  index <- data.frame(sample = c(unique(meta_order$sample),"all"), si = NA, ch = NA, db = NA, ar = NA, fm = NA, nmi = NA)
  meta_order$sample_id <- match(meta_order$sample, unique(meta_order$sample))

  # silhouette neibu
  si <- silhouette(x = meta_order$sample_id, dist = dist(t(matrix_order), method = "manhattan", diag = FALSE, upper = FALSE))
  tem <- data.frame(sample = meta_order$sample, si = si[,3]) %>%
    group_by(sample) %>%
    summarise(si = mean(si))
  index$si <- tem$si[match(index$sample, tem$sample)]

  # Calinski-Harabasz index neibu
  #library(fpc)
  index$ch[nrow(index)] <- calinhara(t(matrix_order), meta_order$sample_id)

  # Davies-Bouldin neibu
  #library(clusterSim)
  index$db[nrow(index)] <- index.DB(t(matrix_order), meta_order$cluster_tem, d=NULL, centrotypes="centroids", p=1, q=2)$DB

  # adjustedRandIndex
  #library(mclust)
  index$ar[nrow(index)] <- adjustedRandIndex(meta_order$sample_id, meta_order$cluster_tem)

  # FM_index
  #library(dendextend)
  index$fm[nrow(index)] <- as.numeric(FM_index(as.character(meta_order$sample_id),as.character(meta_order$cluster_tem)))

  # NMI
  #library(aricode)
  index$nmi[nrow(index)] <- NMI(meta_order$sample_id, meta_order$cluster_tem)

  return(index)
}
