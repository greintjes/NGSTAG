#' NMDS Analysis
#'
#' This function tests whether there is a satistically significant difference between two or more groups of sampling units (community and contextual data)
#' - The function uses normalised community matrix, then does a dissimilartiy analysis using VEGAN vegdist() using the Jaccard index and a subsequent metaMDS() analysis.
#' - It subsequently uses the VEGAN envfit() analysis on a contextual data matrix.
#' - Finally significance test can be done using anosim() and mantel() functions.
#' 
#' @param normalised community matrix, default plot=FALSE
#' @keywords NMDS analysis
#' @export 
#' @examples
#' NGS_Cluster_analysis(x)
#' 
NGS_Significance_function= function(x,y){
  
 
  bray_all<- vegdist(x, method="jaccard")  # calculate Bray-Curtis distance among samples
  mds_bray<<-metaMDS(x, distance="jaccard") # robustness
  ef<= envfit(mds_bray, y, permu=999) 
  
  anosim(x, contex_column_name_1, perm=1000, distance="bray")
  anosim(x, contex_column_name_2, perm=1000, distance="bray")
  anosim(x, contex_column_name_3, perm=1000, distance="bray")
  anosim(x, contex_column_name_4, perm=1000, distance="jaccard")
  anosim(x, contex_column_name_5, perm=1000, distance="jaccard")
 
  comm_dist=vegdist(x, method="bray")
  contex_dist=vegdist(y, method="euc")
  mantel(comm_dist, contex_dist, perm=1000) #  Mantel's permutation test for similarity of two matrices
  # mantel(comm_dist, euk_general_dist, perm=1000) # also works for two communtiy matrixes
  
  
  
  
}

