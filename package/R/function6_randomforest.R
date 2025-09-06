# functions
# randomforest


#' @title Run Random Forest classifier
#' @description Train RF to classify celltypes from VAF matrix and evaluate.
#' @param matrix VAF matrix (variants × cells).
#' @param train_ratio Training proportion.
#' @param meta Metadata with celltype.
#' @param output_path Folder for plots.
#' @param pattern Color palette name.
#' @param celltype_color Named color vector or NULL.
#' @return List with CV metrics and feature importance.
#' @export
run_randomforest <- function(matrix = matrix, train_ratio = 0.7, meta = meta, output_path = getwd(), pattern = "Set1", celltype_color = NULL){
  # check cell
  if( length(unique(meta$cell))==length(unique(meta$cell, colnames(matrix))) ){
    meta <- meta %>%
      filter(cell %in% colnames(matrix)) %>%
      group_by(celltype, cluster_tem)
    meta <- as.data.frame(meta)
    matrix <- t(matrix[,meta$cell])
  }else{
    print("Error: matrix do NOT matches meta\n")
    break
  }

  # make train and test data
  data <- as.data.frame(cbind(as.character(meta$celltype), apply(matrix,2,as.numeric)))
  y <- factor(data[,1])
  X <- as.matrix(apply(data[,-1],2,function(x) as.numeric(x)))
  rownames(X) <- meta$cell
  set.seed(2024)
  index <- caret::createDataPartition(y, p = train_ratio, list = FALSE)
  X_train <- X[index, ]
  X_test <- X[-index, ]
  y_train <- y[index]
  y_test <- y[-index]

  # hyperparameter
  cat("RandomForest: hyperparameters are being trained...\n")
  cv_model_mtry <- get_mtry(X_train = X_train, y_train = y_train) # training mtry
  mtry <- cv_model_mtry$finalModel$mtry
  cv_model_ntree <- get_ntree(X_train = X_train, y_train = y_train, mtry = mtry) # training ntree
  ntree <- get_ntree_best(cv_model_ntree = cv_model_ntree)
  cat(paste0("mtry = ",mtry,"; ntree = ",ntree,"\n"))

  # randomforest
  final_model <- randomForest(x = X_train, y = y_train,
                              ntree = ntree,
                              mtry = mtry,
                              importance = TRUE)
  OOB_err_rate <- signif(min(final_model$err.rate[ntree,1]),3)*100
  cat(paste0("Final OOB estimate of error rate: ", OOB_err_rate, "%\n"))

  # predict
  y_pred <- predict(final_model, X_test)

  # cross validation
  cat("Cross validation: N of significant mutations...\n")
  cv_result = replicate(5, rfcv(data[,-1], factor(as.character(data[,1])), cv.fold = 5), simplify = FALSE)
  cv_error = sapply(cv_result,"[[","error.cv")


  ################################# save figures ###############################
  if(is.null(celltype_color)){
    celltype_color <- get_color(meta = meta, col = "celltype", pattern = pattern)
  }

  names = as.character(unique(meta$celltype))
  dat <- cv_model_mtry$results[,1:2]
  p_mtry <- show_mtry(dat = dat)

  dat <- as.data.frame(cv_model_ntree[[10]]$finalModel$err.rate[,-1])
  dat[is.na(dat)] <- 1
  dat$ntree <- c(1:nrow(dat))
  dat <- reshape2::melt(dat, id.vars = c("ntree"))
  colnames(dat) <- c("ntree", "celltype", "CV_Error")
  p_ntree <- show_ntree(dat = dat, celltype_color = celltype_color)

  dat <- data.frame(celltype = names(randomForest::margin(final_model, y_train)),
                    Prob = as.numeric(randomForest::margin(final_model, y_train))) %>%
    arrange(Prob) %>%
    mutate(index = c(1:length(y_train)))
  p_prob <- show_prob(dat = dat, celltype_color = celltype_color)

  importance <- as.data.frame(importance(final_model))
  dat <- importance %>%
    mutate(Var = rownames(final_model$importance)) %>%
    mutate(Var = ifelse(nchar(Var) > 25, paste0(substring(Var,1,25),".."), Var))
  p_imp_a <- show_imp(dat = dat[,c("Var","MeanDecreaseAccuracy")])
  p_imp_g <- show_imp(dat = dat[,c("Var","MeanDecreaseGini")])

  dat <- cbind(data.frame(n_var = cv_result[[1]]$n.var),cbind(rowMeans(cv_error),cv_error))
  dat <- reshape2::melt(dat, id.vars = c("n_var"))
  colnames(dat) <- c("n_var", "n_cv", "CV_Error")
  p_var <- show_var(dat)
  cat("DONE!\n")

  ################################## output ####################################
  p <- patchwork::wrap_plots(p_mtry, p_ntree, p_prob, p_var, p_imp_a, p_imp_g, byrow = T, ncol = 3)
  ggplot2::ggsave(paste0(output_path,"/","RF_hyperparameter.pdf"), p, width = 23, height = 12)

  importance <- importance[order(importance$MeanDecreaseAccuracy, decreasing = T),]
  confusionMatrix <- confusionMatrix(as.factor(y_pred), as.factor(y_test), mode = "everything")

  list <- list(confusionMatrix = confusionMatrix, importance = importance)
  return(list)
}


#' @title Tune RF mtry
#' @description Cross-validate mtry for Random Forest.
#' @param X_train Feature matrix.
#' @param y_train Response factor.
#' @param n_cv Number of folds.
#' @return caret train object.
get_mtry <- function(X_train = X_train, y_train = y_train, n_cv = 5){
  ctrl <- trainControl(method = "cv", number = n_cv) # Default:5-fold
  grid <- expand.grid(mtry = c(1:(2*ceiling(sqrt(nrow(X_train)))) ))
  cv_model <- train(x = X_train, y = y_train,
                    method = "rf",
                    trControl = ctrl,
                    tuneGrid = grid)
  return(cv_model)
}


#' @title Tune RF ntree
#' @description Cross-validate ntree list for fixed mtry.
#' @param X_train Feature matrix.
#' @param y_train Response factor.
#' @param n_cv Folds.
#' @param mtry Fixed mtry.
#' @param ntree_list Sequence of ntree values.
#' @return List of caret models.
get_ntree <- function(X_train = X_train, y_train = y_train, n_cv = 5, mtry = mtry, ntree_list = seq(100,1000,100)){
  ctrl <- trainControl(method = "cv", number = n_cv) #Default:5-fold
  grid <- expand.grid(mtry = mtry)
  cv_model_ntree <- list()
  for(ntree in ntree_list){
    set.seed(2024)
    cv_model <- train(x = X_train, y = y_train,
                      method = "rf",
                      metric = "Accuracy",
                      trControl = ctrl,
                      tuneGrid = grid,
                      ntree = ntree)
    key <- toString(ntree)
    cv_model_ntree[[key]] <- cv_model
  }
  return(cv_model_ntree)
}


#' @title Select best ntree
#' @description Pick ntree with max accuracy/kappa from CV results.
#' @param cv_model_ntree List from get_ntree.
#' @return Optimal ntree integer.
get_ntree_best <- function(cv_model_ntree = cv_model_ntree){
  results <- resamples(cv_model_ntree)
  accuracy <- summary(results)$statistics$Accuracy
  kappa <- summary(results)$statistics$Kappa
  idx1 <- min(which(accuracy[,"Mean"]==max(accuracy[,"Mean"])))
  idx2 <- min(which(kappa[,"Mean"]==max(kappa[,"Mean"])))
  ntree <- as.numeric(names(cv_model_ntree)[max(idx1,idx2)])
  return(ntree)
}


#' @title Plot mtry curve
#' @description ggplot of CV accuracy vs mtry.
#' @param dat Data frame from caret results.
#' @return ggplot object.
show_mtry <- function(dat = dat){
  p_mtry <- ggplot(dat, aes(x = mtry, y = Accuracy))+
    geom_point()+
    geom_line()+
    theme_bw(base_size = 16)+
    ggtitle("Hyperparameters training: accuracy with mtry")+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_mtry)
}


#' @title Plot ntree curve
#' @description ggplot of CV error vs ntree by celltype.
#' @param dat Long-format data frame.
#' @param celltype_color Named color vector.
#' @return ggplot object.
show_ntree <- function(dat = dat, celltype_color = celltype_color){
  p_ntree <- ggplot(dat, aes(x = ntree, y = CV_Error, color = celltype))+
    geom_point()+
    geom_line()+
    scale_color_manual(values = celltype_color)+
    scale_x_continuous(trans = "log10")+
    theme_bw(base_size = 16)+
    ggtitle("Hyperparameters training: error with ntree")+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_ntree)
}


#' @title Plot margin probability
#' @description ggplot of RF margin vs index.
#' @param dat Data frame with Prob and celltype.
#' @param celltype_color Named color vector.
#' @return ggplot object.
show_prob <- function(dat = dat, celltype_color = celltype_color){
  dat$Prob[which(dat$Prob<0)] <- 0
  p_prob <- ggplot(dat, aes(x = index, y = Prob, color = celltype))+
    geom_point()+
    scale_color_manual(values = celltype_color)+
    theme_bw(base_size = 16)+
    ggtitle("Probability that an observation is judged to be correct")+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_prob)
}


#' @title Plot feature importance
#' @description ggplot of top features colored by mutation type.
#' @param dat Data frame of importance scores.
#' @param top Number of top features.
#' @return ggplot object.
show_imp <- function(dat = dat, top = 50){
  dat <- head(dat[order(dat[,2], decreasing = T),], top)
  mut <- separate(dat, Var, into = c("chr","pos","refalt"), sep = "_")
  mut <- separate(mut, refalt, into = c("ref","alt"), sep = ">")
  mut$type <- paste0(mut$ref,">",mut$alt)
  mut$type[which(mut$type=="G>T")] <- "C>A"
  mut$type[which(mut$type=="G>C")] <- "C>G"
  mut$type[which(mut$type=="G>A")] <- "C>T"
  mut$type[which(mut$type=="A>T")] <- "T>A"
  mut$type[which(mut$type=="A>G")] <- "T>C"
  mut$type[which(mut$type=="A>C")] <- "T>G"
  mut$type[-which(mut$type%in%c("C>A","C>G","C>T","T>A","T>C","T>G"))] <- "others"
  dat <- dat %>% mutate(type =  mut$type[match(dat$Var, rownames(mut))])
  type_color <- c("C>A"="dodgerblue","C>G"="black","C>T"="red3","T>A"="gray50","T>C"="yellowgreen","T>G"="lightpink", "others" = "gold3")

  p_imp <- ggplot(data = dat, aes(x = dat[,2], y = reorder(dat[,1], -dat[,2]), color = type))+
    geom_point(stat='identity', size = 4)+
    theme_bw()+
    scale_color_manual(values = type_color)+
    xlab(colnames(dat)[2])+
    theme(axis.text.y = element_text(size = 16),axis.text.x = element_text(size = 7, angle = 90),axis.title.x = element_blank())+
    coord_flip()
  return(p_imp)
}


#' @title Plot CV error vs variable count
#' @description ggplot of cross-validation error against number of variables.
#' @param dat Long-format data frame with n_var, CV_Error, n_cv.
#' @param n_cv Number of CV folds (for palette).
#' @return ggplot object.
show_var <- function(dat = dat, n_cv = 5){
  getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
  cv_color = getPalette(8)[1:n_cv]
  names(cv_color) <- c(1:n_cv)
  p_var <- ggplot(dat, aes(x = n_var, y = CV_Error, color = n_cv))+
    geom_point()+
    geom_line()+
    theme_bw(base_size = 16)+
    scale_color_manual(values = cv_color)+
    ggtitle("CV: error with number of variables")+
    scale_x_continuous(trans = "log2", breaks = dat$n_var, labels = dat$n_var)+
    theme(plot.title = element_text(hjust = 0.5))
  return(p_var)
}
