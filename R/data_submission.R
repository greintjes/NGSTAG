#' Data submission function
#'
#' This function can be used for data submission of file type .csv. It reads a community data matrix (site by species), a contextual data matrix and a biome (user defined 
#' factors) matrix. The submission uses file.choose() and uses the headers and row.names of the .csv file. 
#' After submission of the data it give total read abundace per sample. Additionally it makes a scatterplot of the contextual data based on a log transformation 
#' 
#' @param input .csv community and environmental matrix
#' @keywords comm, contex, biome
#' @export 
#' @examples
#' NGS_data_submission()
#' 
NGS_data_submission=function(print=TRUE) {
  community <<- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1) # community data file
  class(community)     # type of matrix
  dim(community)       # dimension
  rownames(community)  # row names
  total<<-apply(community, 1, sum)   # total abundance in each sample, should vary due to read abundance variation
  
  # Environmental data submission 
  contextual<<- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1) 
  pairs(contextual)          # plot the data (scatterplot)
  contextual_log <- log10(contextual)     # if some variables look skewed - log transform all variables
  pairs(contextual_log)          # plot the transformed data
  
  png(file="contextual_scatterplot.png", width=1500, height=1450, res=200) # make image of contextual data scatterplot
  pairs(contextual)
  dev.off()
  
  # Factors data submission 
  biome <<- read.csv(file.choose(), sep=";" ,header = TRUE, row.names = 1)
  
}