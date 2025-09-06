#' @title Assign Colors to Categorical Levels
#' @decription Return a named vector of colors for unique values in a metadata column.
#' @param meta Metadata data frame.
#' @param col Column name to color.
#' @param pattern RColorBrewer palette name.
#' @return colors
#' @export
get_color <- function(meta = meta, col = "celltype", pattern = "Set1"){
  id <- which(colnames(meta)==col)
  need_color = length(unique(c(meta[,id])))
  max_color <- brewer.pal.info[pattern,1]
  if(need_color<=max_color){
    getPalette = colorRampPalette(brewer.pal(max_color, pattern))
    need_color = getPalette(max_color)[1:need_color]
  }else{
    getPalette = colorRampPalette(brewer.pal(max_color, pattern))
    need_color = getPalette(need_color)
  }
  names(need_color) <- unique(meta[,id])
  return(need_color)
}


#' @title Color Sub-clusters by Cell Type
#' @decription Generate colors for cluster_new levels, shaded by parental celltype.
#' @param meta Metadata.
#' @param celltype_color Named color vector of celltypes.
#' @param r Blend ratio (0–1) between celltype color and rainbow.
#' @return colors
#' @export
get_color_cluster_new <- function(meta = meta, celltype_color = celltype_color, r = 0.2){
  cluster_new_color <- c()
  for(i in 1:length(celltype_color)){
    j <- celltype_color[i]
    cluster_new_tem <- unique(meta$cluster_new[which(meta$celltype==names(j))])
    if(length(cluster_new_tem)>1){
      num <- length(cluster_new_tem)
      dark <- darken_color(j, darker = -0.1)
      rgb1 <- col2rgb(dark)
      rainbow <- rainbow(num)
      for(i in 1:num){
        rgb2 <- col2rgb(rainbow[i])
        rainbow[i] <- rgb((r*rgb2[1]+(1-r)*rgb1[1])/256,(r*rgb2[2]+(1-r)*rgb1[2])/256,(r*rgb2[3]+(1-r)*rgb1[3])/256)
      }
      cluster_new_color_tem <- rainbow
      names(cluster_new_color_tem) <- cluster_new_tem
      cluster_new_color <- c(cluster_new_color,cluster_new_color_tem)
    }else if(length(cluster_new_tem)==1){
      dark <- darken_color(j, darker = -0.2)
      cluster_new_color_tem <- dark
      names(cluster_new_color_tem) <- cluster_new_tem
      cluster_new_color <- c(cluster_new_color,cluster_new_color_tem)
    }
  }
  return(cluster_new_color)
}


#' @title Darken or Lighten a Color
#' @decription Shift a color toward black or white.
#' @param color Input color.
#' @param darker Positive darkens, negative lightens.
#' @return colors
#' @export
darken_color <- function(color, darker = 0.1) {
  rgb1 <- col2rgb(color)
  if(darker>0){
    rgb2 <- col2rgb("black")
  }else{
    rgb2 <- col2rgb("white")
  }
  k <- abs(darker)
  dark_color <- rgb((k*rgb2[1]+(1-k)*rgb1[1])/256,(k*rgb2[2]+(1-k)*rgb1[2])/256,(k*rgb2[3]+(1-k)*rgb1[3])/256)
  return(dark_color)
}
