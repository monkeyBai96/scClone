# functions for scClone
# bsh2024


#' @title Split FORMAT-like column
#' @description Split colon-separated FORMAT field into named numeric columns.
#' @param dat Input data frame.
#' @param split Column index of values to split.
#' @param name Column index of field names.
#' @param sep Separator character.
#' @return Data frame with named numeric columns.
#' @export
split_FORMAT_col <- function(dat = dat, split = 2, name = 1, sep = ":"){
  cat("Splitting... \n")
  options(warn = -1)
  colnames(dat)[c(name,split)] <- c("name_xxx","split_xxx")
  length = max(str_count(dat[,split], pattern = sep))+1
  split_dat = separate(dat[,c(name,split)], split_xxx, into = paste0("New_",c(1:length)), sep = sep)
  name_dat = separate(dat[,c(split,name)], name_xxx, into = paste0("New_",c(1:length)), sep = sep)
  name_all <- do.call(c,lapply(name_dat[,-1], unique))
  name_all <- as.character(unique(name_all[!is.na(name_all)]))
  pb <- txtProgressBar(min = 0, max = length(name_all), style = 3)
  output <- matrix(NA, ncol = length(name_all), nrow = nrow(dat))
  for(i in 1:ncol(output)){
    setTxtProgressBar(pb, i)
    name_tem <- name_all[i]
    for(j in 2:ncol(name_dat)){
      id <- which(name_dat[,j]==name_tem)
      if(length(id)>0){
        output[id,i] <- split_dat[id,j]
      }
    }
  }
  close(pb)
  cat("DONE! \n")
  output <- as.data.frame(output)
  colnames(output) <- name_all
  return(output)
}


#' @title Split INFO-like column
#' @description Split semicolon-separated INFO field into key=value columns.
#' @param dat Input data frame.
#' @param split Column index to split.
#' @param sep Separator between fields.
#' @param link Separator between key and value.
#' @return Data frame with named columns.
#' @export
split_INFO_col <- function(dat = dat, split = split, sep = ";", link = "="){
  cat("Splitting... \n")
  options(warn = -1)
  dat$name_xxx <- "demo"
  colnames(dat)[c(split)] <- "split_xxx"
  length = max(str_count(dat[,split], pattern = sep))+1
  split_dat = separate(dat[,c(ncol(dat),split)], split_xxx, into = paste0("New_",c(1:length)), sep = sep)
  name_dat <- do.call(cbind,lapply(split_dat, function(x) ifelse(grepl(link,x),sub(paste0(link,".*$"),"",x), x)))
  split_dat <- do.call(cbind,lapply(split_dat, function(x) ifelse(grepl(link,x),sub(paste0("^.*",link),"",x), x)))
  name_all <- do.call(c,lapply(name_dat[,-1], unique))
  name_all <- as.character(unique(name_all[!is.na(name_all)]))
  pb <- txtProgressBar(min = 0, max = length(name_all), style = 3)
  output <- matrix(NA, ncol = length(name_all), nrow = nrow(dat))
  for(i in 1:ncol(output)){
    setTxtProgressBar(pb, i)
    name_tem <- name_all[i]
    for(j in 2:ncol(name_dat)){
      id <- which(name_dat[,j]==name_tem)
      if(length(id)>0){
        output[id,i] <- split_dat[id,j]
      }
    }
  }
  close(pb)
  cat("DONE! \n")
  output <- as.data.frame(output)
  colnames(output) <- name_all
  return(output)
}


#' @title Split VAF column
#' @description Split comma-separated VAF into multiple numeric columns.
#' @param dat Input data frame.
#' @param col Column index to split.
#' @param into New column names vector or NA for auto naming.
#' @param sep Separator character.
#' @return Data frame with split columns.
#' @export
split_VAF_col <- function(dat = dat, col = col, into = NA, sep = ","){
  length = max(str_count(dat[,col], pattern = sep))+1
  if(length(into)==1){
    suppressWarnings({dat <- separate(dat, col, into = paste0(col,"_",c(1:length)), sep = sep)})
  }else{
    if(length(into)==length){
      suppressWarnings({dat <- separate(dat, col, into = into, sep = sep)})
    }else{
      print("Error 'into' provided.\n")
      break
    }
  }
  return(dat)
}


#' @title Change column class
#' @description Convert specified columns to factor or numeric.
#' @param dat Input data frame.
#' @param factor_col Vector of factor column names.
#' @param number_col Vector of numeric column names.
#' @return Data frame with updated classes.
#' @export
change_class <- function(dat = dat, factor_col = factor_col,
                         number_col = setdiff(colnames(dat), factor_col)){
  dat <- as.data.frame(dat)
  for(i in 1:ncol(dat)){
    col_tem <- colnames(dat)[i]
    if(col_tem %in% factor_col){
      dat[,i] <- factor(dat[,i])
    }else if(col_tem %in% number_col){
      dat[,i] <- as.numeric(dat[,i])
    }else{
      break
      cat("Error 'col' provided.\n")
    }
  }
  return(dat)
}


#' @title First column to rownames
#' @description Move first column to rownames and return data frame.
#' @param data Input data frame.
#' @return Data frame with updated rownames or error if duplicates.
#' @export
col1_to_rowname <- function(data = data){
  data <- as.data.frame(data)
  rownames <- data[,1]
  if(any(duplicated(rownames))==T){
    cat("ERROR: Duplicated rownames provided.\n")
  }else{
    dat <- as.data.frame(data[,-1])
    colnames(dat) <- colnames(data)[-1]
    rownames(dat) <- rownames
    return(dat)
  }
}


#' @title Read ANNOVAR outputs
#' @description Batch read ANNOVAR txt files and merge into one table.
#' @param path Folder containing *.txt files.
#' @param col_keep Column indices to retain.
#' @param col_name Optional new column names vector.
#' @return Merged data frame with cell column.
#' @export
read_annovar <- function(path = path, col_keep = col_keep, col_name = NULL){
  files <- list.files(path = path, pattern = ".txt")
  file_names <- as.character(substring(files,1,nchar(files)-19))
  tag <- as.character(substring(files,nchar(files)-18,nchar(files)))[1]
  cat(paste0("Reading ", length(file_names)," files","...\n"))
  pb <- txtProgressBar(min = 0, max = length(file_names), style = 3)
  file_list <- list()
  for(i in 1:length(file_names)){
    raw <- read.table(paste0(path,"/",file_names[i], tag), header = T, quote = "", sep = "\t")
    raw <- raw[,col_keep]
    raw$cell <- file_names[i]
    file_list[[i]] <- raw
    setTxtProgressBar(pb, i)
  }
  data <- do.call(rbind, file_list)
  if(!is.null(col_name)){
    colnames(data) <- col_name
  }
  close(pb)
  cat("DONE!\n")
  return(data)
}
