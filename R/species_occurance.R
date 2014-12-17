#' Species_occurance_function
#'
#' This function anaylsis a normalised community matrix. 
#' - It takes a species by site matrix and analyses species occurance within and between sites. 
#' - Optional transformation of the data matrix to emphasise or deemphasis dominants 
#' - How many sites does each species occur in ? This part of the function first calls (veg=0) 
#'   which evaluates True/flase or 1/0, and then sums the column values
#' - Species occurance all sites or 1 or less sites (spc.pres X, X being number of sites) 
#' - Gets the average cover for each species (sum in site divided by the number of sites) 
#' - optional identify which are the dominant organism using listpts function
#' - Calculate species area relationship 
#' - all results are stored in file speceis_occurance_results.csv includes names of species in 
#'   all sites or in 1 or less sites. Also top 10 abundance species are included. All graphs 
#'   (.png images) are stored in the working directory. They are high quailty images 
#' 
#' @param normalised community matrix, echo=TRUE
#' @keywords species occurance
#' @export 
#' @examples
#' NGS_species_occurance_functiom(x, echo=TRUE)
#' 
NGS_species_occurance_function<-function(x, echo=T){
  
  
  # How many sites does each species occur in.
  spc.pres<-apply(x>0,2,sum)    # to get number of presences for each species. (1,0 presence and absence)
  png(filename="Species_Occurance.png", width=1500, height=1450, res=250)
  plot(sort(spc.pres), main="Species Occurance", xlab='Cumulative Number of Genera',ylab='Number of Sites', col="black", pch=16)
  dev.off()
  png(filename="Species_Occurance(log).png", width=1500, height=1450, res=250)
  plot(sort(spc.pres),log='x',main="Species Occurance(log)", xlab='Cumulative Number of Genera (log)',ylab='Number of Sites', col="black", pch=16)   				#if skewed put y axis on log scale
  dev.off()
  # Abundant species (enter number of site (e.g. 35 or more))  
  spc.pres29<-	spc.pres[spc.pres>=29]
  # Rare species (that occur 1 or less sites)
  spc.pres1<-	spc.pres[spc.pres<=1]    
  
  # mean abundance of each species  
  tmp <- 	apply(x,2,sum)	
  spc.mean <- tmp/spc.pres     
  png(filename="Mean Coverage When a Genera Occurs.png", width=1500, height=1450, res=250)
  plot(sort(spc.mean),main="Mean Coverage When a Genera Occurs", xlab="Cumulative Number of Sites",ylab="Mean Abundance", col="black", pch=16)
  dev.off()
  # mean abundance of each species against number of sites it occur in
  png(filename="Species Abundance (high).png", width=1500, height=1450, res=250)
  plot(spc.pres,spc.mean, main="Cumulative Count of Genera Against Mean Abundance", xlab="Cumulative Count of Sites", ylab="Mean Abundance", col="black", cex=1, pch=16)	
  dev.off()
  
  # optional identify which are the dominant organism using listpts function
  # need to click near points of interest and use ESC to stop identification of points
  plot(spc.pres,spc.mean, main="Cumulative Count of Genera Against Mean Abundance", xlab="Cumulative Count of Sites", ylab="Mean Abundance", col="grey60", cex=1)
  listpts <- identify(spc.pres,spc.mean,names(x), cex=0.3) 	# list identify abundance species
  dev.copy(png, filename="Species Abundance with Names (high).png", width=1500, height=1450, res=250)
  dev.off()
  
  plot(spc.pres,spc.mean, main="Cumulative Count of Genera Against Mean Abundance", xlab="Cumulative Count of Sites", ylab="Mean Abundance", col="grey60", cex=1)
  listpts2 <- identify(spc.pres,spc.mean,names(x), cex=0.3) 	# list identify abundance species
  dev.copy(png, filename="Species Abundance with Names(low).png", width=1500, height=1450, res=250)
  dev.off()  		
  
  # Calculate species area relationship
  plt.pres<-apply(x>0,1,sum) # to calculate the number of species in each site (for 1 means for column)
  plt.sum <- apply(x,1,sum) 	# to calculate the total number of species on each site
  png(filename="Relationship between number of Genera per site and total area.png", width=1500, height=1450, res=250)
  plot(plt.pres,plt.sum, main="Number of Genera  Against Abundance of Genera ", xlab="Number of Genera per Site", ylab="Abundance of Genera per Site", col="black")  	# see the relationship between number of species/plot and total number of species 
  res=lm(plt.sum~plt.pres)
  res
  abline(res) #draw line for relationship between species/area
  dev.off()
  # calculate top 10 species from occurance  
  ctmp<-apply(x,2, sum)
  ctmp<-sort(ctmp, decreasing=TRUE)
  top10<-ctmp[1:10]
  
  results<-c(spc.pres29=spc.pres29, spc.pres1=spc.pres1, listpts=listpts, top10=top10)
  write.table(results, file="speceis_occurance_results.csv", sep = "\t", row.names=TRUE, col.names=TRUE)
  return(results)
}
