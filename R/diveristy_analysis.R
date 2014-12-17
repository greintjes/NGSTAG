#' Diveristy analysis of NGS data
#'
#' This function anaylsis a normalised community matrix using Vegan based functions. It requires a species by site matrix and analyses this for species RICHNESS, EVENESS and DIVERSITY INDECIES ANALYSIS
#' - Does species analysis :richness 
#' - Does species frequency analysis (freuqency of occurance)
#' - Does Shannon Weaver analysis 
#' - Does Simpson analysis 
#' - Does Inverse Simpsoms analysis 
#' - Does Evenness analysis 
#' - All results are stored in file diversity_results.csv, all graphs (.png images) are stored 
#'   in the working directory. They are high quailty images. 
#' 
#' @param normalised community matrix, print default =TRUE
#' @keywords diveristy analysis
#' @export 
#' @examples
#' NGS_diversity_index_function(x, print=TRUE)
#' 

NGS_diversity_indecies_function<-function(x, print=TRUE)
{
  
  #number of species
  sn<-specnumber(x)   					
  png(filename="Richness.png", width=1500, height=1450, res=250)
  par( bty="n")		
  plot(sn, main="Richness", col="black", xlab="", ylab="") # no x and y labels, x and y limits set.
  #lines(sn,col="grey20",lwd=3) 		        		# join points with pint line, width 3
  mtext(side=1,line=2.5,cex=1,"Station") 			# add text under side 1 of plot
  mtext(side=2,line=2.5,cex=1,"Richness") 		# add text under side 2 of plot	
  text(sn, row.names(x), cex=0.6, pos=3, col="grey20") 	# add text to dots, pos is location of text
  dev.off()
  
  #frequencies of species
  snfreq<- specnumber(x, MARGIN=2) 			
  png(file="Genus_Frequency.png", width=1500, height=1450, res=250)
  par( bty="n")					
  plot(snfreq, main="Genus Frequency", col="black", xlab="", ylab="")
  mtext(side=1,line=2.5,cex=1,"Genus Frequencey") 
  mtext(side=2,line=2.5,cex=1,"Station") 			      
  dev.off()
  
  # shannon-weaver
  h<-diversity(x)						
  png(file="Shannon_Weaver.png", width=1500, height=1450, res=250)	
  par( bty="n")						
  plot(h, main="Shannon-Weaver", col="black", xlab="", ylab="") 
  mtext(side=1,line=2.5,cex=1,"Station") 			    
  mtext(side=2,line=2.5,cex=1,"Shannon Index") 	
  text(h, row.names(x), cex=0.6, pos=1, col="grey20")
  dev.off()
  
  # simpsons
  s<-diversity(x,index="simpson")
  png(file="Simpson.png", width=1500, height=1450, res=250)
  par( bty="n")	
  plot(s, main="Simpsons", col="black", xlab="", ylab="")
  mtext(side=1,line=2.5,cex=1,"Station") 			
  mtext(side=2,line=2.5,cex=1,"Simpsons Index") 		
  text(s, row.names(x), cex=0.6, pos=1, col="grey20")
  dev.off()
  
  # inverse simpsons 
  i<-diversity(x, index="invsimpson")
  png(file="Inverse_Simpson.png", width=1500, height=1450, res=250)
  par( bty="n")	
  plot(i, main="inverse Simpsons", col="black", xlab="", ylab="")
  mtext(side=1,line=2.5,cex=1,"Station") 			
  mtext(side=2,line=2.5,cex=1,"Inverse Simpsons Index") 
  text(i, row.names(x), cex=0.6, pos=1, col="grey20")
  dev.off()
  
  #Evenness
  j<-h/log(specnumber(x))
  png(file="Eveness.png", width=1500, height=1450, res=250)
  par( bty="n")	
  plot(j, main="Evenness", col="black" , xlab="", ylab="") # ylim may cause issue change to just plot(j)
  mtext(side=1,line=2.5,cex=1,"Station") 			
  mtext(side=2,line=2.5,cex=1,"Eveness") 			
  text(j, row.names(x), cex=0.6, pos=1, col="grey20")
  dev.off()

  
  #Results return values for all diversity indecies in diversity_results.csv file
  results<-list(sn=sn, h=h, s=s, i=i, j=j)
  write.table(results, file = "diversity_results.csv", sep = "\t", row.names=TRUE, col.names=TRUE)
  return(results)
  
}
