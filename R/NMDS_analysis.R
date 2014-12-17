#' NMDS Analysis
#'
#' This function anaylsis a normalised community matrix, contextual matrix and factor matrix
#' - It requires a species by site matrix, contextual data matrix and selection factors matrix
#' - It calculates a distance matrix using vegdist() by bray curtis and jaccard dissimilarity analysis.
#' - Then it uses the metaMDS() function to perform nonmetric multidimensional scaling (NMDS).
#' - Subsequently it prints a stress plot to assess goodness of the ordination fit 
#' - Then it plots a basic NMDS, this can be done with "sites" or "species" labeled.
#' - Subsequently more complex NMDS plots can be made with seperations or coloring based on user defined factors or contextual data.
#' - Then the funtion does rankindex() of community and contex data and analysis the significant of environmental factors on the community dissimilarity using the envfit() function.
#' - Results are listed after analysis
#' 
#' @param normalised community matrix, default plot=TRUE, default site=TRUE
#' @keywords NMDS analysis
#' @export 
#' @examples
#' NGS_NMDS_function(x,y,z, plot=FALSE)
NGS_NMDS_function<- function(x,y,z, plot=TRUE, sites=TRUE){
  
  # Calculate Bray-Curtis distances among samples based on species P/A/abundance 
  bray_all <- vegdist(x, method="bray")  
  mds_bray<-metaMDS(x, distance="bray") # robustness
   png(file="stressplot.png", width=1500, height=1450, res=200)   # Assess goodness of ordination fit (stress plot)
   stressplot(mds_bray, bray_all, pch=".", p.col="black", l.col="black") 
   dev.off()
   
  jaccard_all <- vegdist(x, method="jaccard")  
  mds_jaccard<-metaMDS(x, distance="jaccard")
    png(file="stressplot_jaccard.png", width=1500, height=1450, res=200)
    stressplot(mds_jaccard, jaccard_all, pch="." , p.col="black", l.col="black")  
    dev.off()
  
  # Plot NMDS basic plot  
  if (sites ==TRUE){
  png(file="NMDS_Basic_vegan.png", width=1500, height=1450, res=200)
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  plot(mds_bray, display = "sites",   type = "p", xlim=c(-1,1), ylim=c(-1, 1))       # plot site scores as text
  points(mds_bray, display = "sites", cex = 0.8, pch=21, col="black", bg="grey80")
  #text(mds_bray, display = "sites", cex=0.7, col="grey20", pos=2)
  
  dev.off()  
  
  }
  
  if (sites ==FALSE) {
    png(file="NMDS_Basic_vegan.png", width=1500, height=1450, res=200)
    par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
    plot(mds_bray, display = "species",   type = "p", xlim=c(-1,1), ylim=c(-1, 1))       # plot site scores as text
    points(mds_bray, display = "species", cex = 0.8, pch=21, col="black", bg="grey80")
    #text(mds_bray, display = "sites", cex=0.7, col="grey20", pos=2)
    
    dev.off()  }
    
  

  #Plot NMDS showing specific factor.
  
  png(file="NMDS_size_fractions.png", width=1500, height=1450, res=200)
  
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-1,1), ylim=c(-1, 1), cex.axis=1)   # don't plot anything yet
  points(mds.fig, "sites", pch = 19, col = "green", select = factor_1 == 
           "Factor_1_title_1")     # plot just the samples, colour by habitat, pch=19 means plot a circle
  points(mds.fig, "sites", pch = 19, col = "blue", select = factor_1 == 
           "Factor_1_title_2")
  
  ordispider(mds_bray, factor_1, label = TRUE, cex=1.5) # add confidence ellipses around habitat types
  dev.off()
  
  
  # NMDS with meansure (environmental) parameters.
  rank<-rankindex(y, x, c("euc","man","bray","jac","kul")) 
  ef<-envfit(mds_bray, y, permu=999) # environmental parameter significance
  
  if (plot ==TRUE){
    png(file="NMDS_vegan_contex.png", width=1500, height=1450, res=200)
    plot(mds_bray, display="sites",type = "p", xlim=c(-1,1), ylim=c(-1,1))
    points(mds_bray, display = "sites", cex = 0.8, pch=21, col="black", bg="grey80")
    text(mds_bray, display = "sites", cex=0.7, col="grey20", pos=1)
    plot(ef,p.max=0.05,cex=0.6, col="black")
    dev.off()
    
   
  }
  
  results<-list(rank=rank, ef=ef, ef2=ef2)
  return(results)
  
}
