#' Cluster analysis
#'
#' This function perfroms a dissimilarity analysis using the vegdist function (VEGAN) on a normalised community matrix. It subsequently performs a hierarchical cluster analysis via the (stats) hclust function.
#' -The function requires a species by site matrix
#' -It makes labels based on row.names
#' -Trees are calulated using Jaccard or bray curtis similarity matrixes vegdist() 
#'   and clustered by "single", "average" or "complete" linkage algorithms hclust() functions. 
#'   Trees are then plotted together for comparison. 
#' - The function directly plots the anylsis.
#' 
#' @param normalised community matrix, default na.rm=TRUE
#' @keywords cluster analysis
#' @export 
#' @examples
#' NGS_Cluster_analysis(x)
#' 
NGS_Cluster_analysis<-function(x){
  
  
  # bray curtis tree calculation
  bray_cluster<-vegdist(x, "bray", na.rm=TRUE)  # calculate Bray-Curtis distance among samples
  treesin=hclust(bray_cluster, method="single") 
  treeave=hclust(bray_cluster, method="average")  # cluster communities using average-linkage algorithm
  treecom=hclust(bray_cluster, method="complete")
  
  png(file="Tree_braycurtis.png", width=1500, height=1450, res=200)
  par(mfrow=c(3,1))
  plot(treesin, hang=-1, main="single")
  plot(treeave, hang=-1, main="average")  # cluster communities using average-linkage algorithm
  plot(treecom, hang=-1, main="complete")
  dev.off()
  
  # Jaccard tree calculation
  jaccard_cluster<-vegdist(x, "jaccard", na.rm=TRUE)  # calculate Jaccard distance among samples
  treesin=hclust(jaccard_cluster, method="single")
  treeave=hclust(jaccard_cluster, method="average")
  treecom=hclust(jaccard_cluster, method="complete")
  
  png(file="Tree_jaccard.png", width=1500, height=1450, res=200)
  par(mfrow=c(3,1))
  plot(treesin, hang=-1, main="single")
  plot(treeave, hang=-1, main="average")
  plot(treecom, hang=-1, main="complete")
  dev.off()
  
  
}
