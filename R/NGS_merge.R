#' Merging function
#'
#' This function allows you to merge several dataframes based on their taxonomy (colname). This function is based on the rbind.all.columns() from http://amywhiteheadresearch.wordpress.com/2013/05/13/combining-dataframes-when-the-columns-dont-match/. It finds the difference between the column names between two vectors (setdiff()) and then adds dummy columns to each vector representing the missing columns and adds NA into the column as data. Then it uses rbind() to combines the vectors. 
#' @param c_merge_taxonomy(x,y,)  X =species site vector 1, y= species site vector 2. The specific function takes two vectors and combines their taxonomy and make a new vector called merge1. this can be done multiple times to add several vector togehter.
#' @param d_sort.removeNA(x), X= merge vector with NA. NA are replaced with 0's and the vector is ordered.
#' @export 
#' @keywords merge
#' NGSmerge(x,y)

NGSmerge <- function(x,y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  y[, c(as.character(x.diff))] <- NA
  
  merge<<-return(rbind(x, y))
  
  #remove NA s
  merge[is.na(merge)] <- 0.0000000
  # order y species names
  All_sites_merged_ordered<<-merge[,order(names(merge))]
  write.table(All_sites_merged_ordered, file="all_merged_ordered.csv", sep=";")
  
}