# functions for scClone
# bsh2024


#' @title Filter mutations
#' @description Iteratively filter cells and mutations by thresholds.
#' @param mut Mutation data frame.
#' @param min_celltype Min cells per cell type.
#' @param min_mut_per_cell Min mutations per cell.
#' @param min_cell_one_mut Min cells sharing one mutation.
#' @param meta Metadata data frame.
#' @return Filtered mutation data frame.
#' @export
mutation_filter <- function(mut = mut, min_celltype = 10, min_mut_per_cell = 10, min_cell_one_mut = 10, meta = meta){
  single <- mut %>% filter(cell %in% meta$cell)
  n_cell <- length(unique(single$cell))
  n_var <- length(unique(single$var))
  n_mut <- nrow(single)

  # cell filter 1
  table_cell <-  as.data.frame(table(single$cell))
  table_cell$celltype <- meta$celltype[match(table_cell$Var1, meta$cell)]
  table_cell$label <- ifelse(table_cell$Freq>min_mut_per_cell, "keep", "no")
  table_cell <- table_cell %>% arrange(label, celltype, -Freq)
  # keep rare cell
  rare_cell <- as.data.frame(table(table_cell$celltype, table_cell$label)) %>% filter(Var2=="keep")
  for(i in 1:nrow(rare_cell)){
    if(rare_cell$Freq[i]<min_celltype){
      n_more <- min(min_celltype-rare_cell$Freq[i], nrow(table_cell %>% filter(celltype==rare_cell$Var1[i] & label=="no")))
      idx <- head(which(table_cell$celltype==rare_cell$Var1[i] & table_cell$label=="no"),n_more)
      if(length(idx)>0){
        table_cell$label[idx] <- "rare"
      }
    }
  }
  cell_high <- as.character(table_cell$Var1[which(table_cell$label!="no")])
  single <- subset(single, cell%in%cell_high)
  # mutation filter 1
  table_var <- as.data.frame(table(single$var))
  mut_high <- as.character(table_var$Var1[which(table_var$Freq > min_cell_one_mut)])
  single <- single %>% filter(var %in% mut_high)
  # cell filter 2
  table_cell <-  as.data.frame(table(single$cell))
  table_cell$celltype <- meta$celltype[match(table_cell$Var1, meta$cell)]
  table_cell$label <- ifelse(table_cell$Freq>min_mut_per_cell, "keep", "no")
  table_cell <- table_cell %>% arrange(label, celltype, -Freq)
  # keep rare cell
  rare_cell <- as.data.frame(table(table_cell$celltype, table_cell$label)) %>% filter(Var2=="keep")
  for(i in 1:nrow(rare_cell)){
    if(rare_cell$Freq[i]<min_celltype){
      n_more <- min(min_celltype-rare_cell$Freq[i], nrow(table_cell %>% filter(celltype==rare_cell$Var1[i] & label=="no")))
      idx <- head(which(table_cell$celltype==rare_cell$Var1[i] & table_cell$label=="no"),n_more)
      if(length(idx)>0){
        table_cell$label[idx] <- "rare"
      }
    }
  }
  cell_high <- as.character(table_cell$Var1[which(table_cell$label!="no")])
  single <- subset(single, cell%in%cell_high)
  # mutation filter 2
  table_var <- as.data.frame(table(single$var))
  mut_high <- as.character(table_var$Var1[which(table_var$Freq > min_cell_one_mut)])
  single <- single %>% filter(var %in% mut_high)

  n_cell_s <- length(unique(single$cell))
  n_var_s <- length(unique(single$var))
  n_mut_s <- nrow(single)
  cat(paste0(n_mut_s,"/",n_mut, " mut keeped. ",n_cell_s,"/",n_cell, " cell keeped. ",n_var_s,"/",n_var, " var keeped.\n"))
  return(single)
}


#' @title Beta density
#' @description Evaluate beta density at given points.
#' @param x Numeric vector of quantiles.
#' @param a First shape parameter.
#' @param b Second shape parameter.
#' @return Density values.
#' @export
beta_function <- function(x,a,b){
  dbeta(x, shape1 = a, shape2 = b)
}


#' @title Transform allele fraction
#' @description Apply nonlinear transformation to allele fraction.
#' @param x Numeric vector in (0,1).
#' @param AB Transformation parameter.
#' @return Transformed values.
#' @export
transform_function <- function(x, AB){
  #VAF_tem <- x
  ## mapping
  #mapping_function <- function(x){
  #  y <- (x-VAF_tem)*(AB+0.5)/(x-0.5) # F1
  #  return( y + sqrt(0.25-(VAF_tem-0.5)^2) - 0.5) # F2 ÏÂ°ëÔ²
  #}
  ## solve
  #if(VAF_tem!=0.5){
  #  solve <- rootSolve::multiroot(f = mapping_function, start = c(VAF_tem))
  #  return(solve$root)
  #}else{
  #  return(0.5)
  #}
  root <- 0.5 + ((x-0.5)*(AB+0.5))/(AB+sqrt(0.25-(x-0.5)^2))
  return(root)
}


#' @title Infer AB parameter
#' @description Estimate transformation parameter from alt/ref counts and capture rate.
#' @param sum_alt Alt read counts.
#' @param sum_ref Ref read counts.
#' @param capture Capture efficiency.
#' @return AB value.
#' @export
infer_AB <- function(sum_alt,sum_ref,capture){
  a = max(sum_alt,sum_ref)
  b = min(sum_alt,sum_ref)
  res = (a-b)*capture/((b+1)^3)
  res <- exp(-res)
  return(res)
}


#' @title Infer true VAF
#' @description Integrate transformed beta posterior to estimate allele fraction.
#' @param ALT Alt read count.
#' @param REF Ref read count.
#' @param AB Transformation parameter.
#' @param lower Integration lower bound.
#' @param upper Integration upper bound.
#' @return Estimated VAF.
#' @export
infer_vaf <- function(ALT = ALT, REF = REF, AB = AB, lower = 0, upper = 1){
  expectation <- cubature::adaptIntegrate(f = function(x){transform_function(x, AB=AB)*beta_function(x, a=ALT+1, b=REF+1)},lower = lower,upper = upper)$integral
  #second_moment <- cubature::adaptIntegrate(f = function(x){transform_function(x,AB=AB)^2*beta_function(x, a=ALT+1, b=REF+1)},lower = lower,upper = upper)$integral
  #variance <- second_moment - expectation^2
  expectation <- round(expectation, 2)
  #variance <- round(variance, 2)
  return(c(expectation))
}


#' @title Build VAF matrix
#' @description Estimate VAF from REF/ALT count matrices via transformation.
#' @param matrix_REF REF count matrix.
#' @param matrix_ALT ALT count matrix.
#' @param fast_cutoff Proportion threshold for fast 0/1/2 encoding.
#' @return VAF matrix.
#' @export
make_VAF_matrix <- function(matrix_REF, matrix_ALT, fast_cutoff = 0.9){
  if(all(dim(matrix_REF)==dim(matrix_ALT))){
    cat(paste0(dim(matrix_REF)[1]," x ",dim(matrix_REF)[2]," matrix provided.\n"))
  }else{
    stop("ERROR: Different matrix provided.\n")
  }
  pb <- txtProgressBar(min = 0, max = nrow(matrix_ALT), style = 3)
  options(warn = -1)
  matrix_VAF <- matrix(NA, nrow = nrow(matrix_REF), ncol = ncol(matrix_REF))
  for(i in 1:nrow(matrix_VAF)){
    setTxtProgressBar(pb, i)
    alt_tem <- as.numeric(matrix_ALT[i,])
    ref_tem <- as.numeric(matrix_REF[i,])
    pure_alt <- length(which(alt_tem>0)) / length(alt_tem)
    if(pure_alt > fast_cutoff){
      matrix_VAF[i,] <- ifelse(alt_tem>0,ifelse(ref_tem==0, 2, 1), NA) # fast it
    }else{
      capture <- length(which(ref_tem+alt_tem > 0))/length(ref_tem)
      AB <- infer_AB(sum(alt_tem[!is.na(alt_tem)]),
                     sum(ref_tem[!is.na(ref_tem)]),
                     capture)
      for(j in 1:ncol(matrix_VAF)){
        if(as.numeric(matrix_ALT[i,j])+as.numeric(matrix_REF[i,j]>0)){
          matrix_VAF[i,j] <- infer_vaf(ALT = as.numeric(matrix_ALT[i,j]), REF = as.numeric(matrix_REF[i,j]), AB = AB)
        }
      }
    }
  }
  close(pb)
  cat("DONE! \n")
  return(matrix_VAF)
}


#' @title Round VAF matrix
#' @description Discretize VAF to 0/1/2 states.
#' @param matrix VAF numeric matrix.
#' @return Integer matrix.
#' @export
rounding_VAF <- function(matrix) {
  matrix[matrix < 0.2] <- 0
  matrix[matrix > 1.2] <- 2
  matrix[(matrix >= 0.2) & (matrix <= 1.2)] <- 1
  return(matrix)
}
