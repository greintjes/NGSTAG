

NGS_data_submission=function(print=TRUE) {
    comm <<- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1) # community data file
      class(comm)     # type of matrix
      dim(comm)       # dimension
      rownames(comm)  # row names
      apply(comm, 1, sum)   # total abundance in each sample, should vary due to read abundance variation

# Environmental data submission 
  contex<<- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1) 
      pairs(contex)          # plot the data (scatterplot)
      contex_log <- log10(contex)     # if some variables look skewed - log transform all variables
      pairs(contex_log)          # plot the transformed data
      
      png(file="Contex_scatterplot.png", width=1500, height=1450, res=200) # make image of contextual data scatterplot
      pairs(contex)
      dev.off()

# Factors data submission 
  biome <<- read.csv(file.choose(), sep=";" ,header = TRUE, row.names = 1)

  }

 
NGS_normalise_abundance_function =function(x){  
  comm_stand<<- decostand(comm, method = "total")  # "total" divide by sample total abundance, "norm" make margin sum of 
                                             # squares equal to one
  norm_check<<-apply(comm, 1, sum)            # check total abundance in each sample, should be all 1 now.

 

}

#

NGS_diversity_index_function<-function(x, print=TRUE)
{
  
#species number 
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
    png(file="Eveness_bacteria.png", width=1500, height=1450, res=250)
    par( bty="n")	
    plot(j, main="Evenness", col="black" , xlab="", ylab="") # ylim may cause issue change to just plot(j)
    mtext(side=1,line=2.5,cex=1,"Station") 			
    mtext(side=2,line=2.5,cex=1,"Eveness") 			
    text(j, row.names(x), cex=0.6, pos=1, col="grey20")
    dev.off()
#species richness per factor (province)
    attach(biome)# attach dataframe with factor!!
    png(file="Species_richness_between_habitats.png", width=1500, height=1450, res=200)
    boxplot(specnumber(comm) ~ province, xlab = "number of species")
    dev.off()

#Results return values for all diversity indecies in diversity_results.csv file
  results<-list(sn=sn, h=h, s=s, i=i, j=j)
  write.table(results, file = "diversity_results.csv", sep = "\t", row.names=TRUE, col.names=TRUE)
  return(results)
   
}



NGS_species_occurance_functiom<-function(x, echo=T){
  

# How many sites does each species occur in.
  spc.pres<-apply(x>0,2,sum)    # to get number of presences for each species. (1,0 presence and absence)
    png(filename="Species_Occurance.png", width=1500, height=1450, res=250)
    plot(sort(spc.pres), main="Species Occurance", xlab='Cumulative Number of Genera',ylab='Number of Sites', col="black", pch=16)
    dev.off()
    png(filename="Species_Occurance(log).png", width=1500, height=1450, res=250)
    plot(sort(spc.pres),log='x',main="Species Occurance(log)", xlab='Cumulative Number of Genera (log)',ylab='Number of Sites', col="black", pch=16) 					#if skewed put y axis on log scale
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
    png(filename="De Relationship between number of Genera per site and total area.png", width=1500, height=1450, res=250)
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



NGS_Cluster_analysis<-function(x){
  station =row.names(x)
  
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
# bootstrapp tree   
  t_vector=t(x)     # transpose dataset
  tree=pvclust(t_vector, method.hclust="ward", nboot=1000, method.dist="euclidean")
    png(file="Tree_bootstrapped.png", width=1500, height=1450, res=200)
    par(mfrow=c(1,1))
    plot(tree, hang=-1, main="boot") 								
    pvrect(tree, alpha=.95, border=2, lwd=3)		# bootstrap							
    dev.off()
  
}



NGS_NMDS_function<- function(x,y,z, plot=FALSE){
  
# Calculate Bray-Curtis distances among samples based on species P/A/abundance 
  bray_all <- vegdist(x, method="bray")  
  mds_bray<-metaMDS(x, distance="bray") # robustness
    png(file="stressplot.png", width=1500, height=1450, res=200)   # Assess goodness of ordination fit (stress plot)
    stressplot(mds_bray, bray_all, pch=".", p.col="black", l.col="black") 
    dev.off()
# Calculate Bray-Curtis distances among samples based on species P/A/abundance   
#  jaccard_all <- vegdist(x, method="jaccard")  
#  mds_jaccard<-metaMDS(x, distance="jaccard")
#    png(file="stressplot_jaccard.png", width=1500, height=1450, res=200)
#    stressplot(mds_jaccard, jaccard_all, pch="." , p.col="black", l.col="black")  
#    dev.off()
  
# Plot NMDS basic plot  
  png(file="NMDS_Basic_vegan.png", width=1500, height=1450, res=200)
par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
    plot(mds_bray, display = "sites",   type = "p", xlim=c(-1,1), ylim=c(-1, 1))       # plot site scores as text
    points(mds_bray, display = "sites", cex = 0.8, pch=21, col="black", bg="grey80")
    #text(mds_bray, display = "sites", cex=0.7, col="grey20", pos=2)
tmp<-with(contex, ordisurf(mds_bray, Temperature, add=TRUE, col="blue")) 

  dev.off()  

#Plot NMDS showing specific factors ( size fraction with spider)  

  png(file="NMDS_size_fractions.png", width=1500, height=1450, res=200)

par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
    mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-1,1), ylim=c(-1, 1), cex.axis=1)   # don't plot anything yet
    points(mds.fig, "sites", pch = 19, col = "green", select = size_fraction == 
         "Free living")   # plot just the samples, colour by habitat, pch=19 means plot a circle
    points(mds.fig, "sites", pch = 19, col = "blue", select = size_fraction == 
         "Attached")

    ordispider(mds_bray, size_fraction, label = TRUE, cex=1.5) # add confidence ellipses around habitat types
  dev.off()
  
#Plot NMDS showing specific factors (province with legend)  
  # col=c("aquamarine", "chartreuse",  "royalblue", "grey80",   "tomato", "coral", "gold", "limegreen", "pink",
  # "lightskyblue4","mediumturquoise","mediumseagreen")

 # col1=c("aquamarine", "chartreuse",  "royalblue", "grey80",   "tomato", "gold")

#  png(file="NMDS_Province.png", width=1500, height=1450, res=200)
 #   mds.fig <- ordiplot(mds_jaccard, type = "none")   # don't plot anything yet
 #   points(mds.fig, "sites", pch = 19, col = "aquamarine", select = province == 
 #          "NADR")   # plot just the samples, colour by habitat, pch=19 means plot a circle
 #   points(mds.fig, "sites", pch = 19, col = "chartreuse", select = province == 
 #          "NAST")
 #   points(mds.fig, "sites", pch = 19, col = "royalblue", select = province == 
 #          "NATR")
 #   points(mds.fig, "sites", pch = 19, col = "grey80", select = province == 
 #          "WTRA")
 #   points(mds.fig, "sites", pch = 19, col = "tomato", select = province == 
 #          "SATL")
 #   points(mds.fig, "sites", pch = 19, col = "gold", select = province == 
 #          "SSTC")
 #   legend("bottomright", legend = c("NADR", "NAST","NATR","WTRA","SATL","SSTC"), # select cluster number 
 #        col = col1 ,  pt.bg = col1, bty = "n", pch = 21)
 # dev.off()
  
#Plot NMDS showing specific factors (province with spider) 
  png(file="NMDS_Province_spider.png", width=1500, height=1450, res=200)
    mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-1,1), ylim=c(-1, 1))   # don't plot anything yet
    points(mds.fig, "sites", pch = 19, col = "aquamarine", select = province == 
           "NADR")   # plot just the samples, colour by habitat, pch=19 means plot a circle
    points(mds.fig, "sites", pch = 19, col = "chartreuse", select = province == 
           "NAST")
    points(mds.fig, "sites", pch = 19, col = "royalblue", select = province == 
           "NATR")   # plot just the samples, colour by habitat, pch=19 means plot a circle
    points(mds.fig, "sites", pch = 19, col = "grey80", select = province == 
           "WTRA")
    points(mds.fig, "sites", pch = 19, col = "tomato", select = province == 
           "SATL") 
    points(mds.fig, "sites", pch = 19, col = "gold", select = province == 
           "SSTC")
    ordispider(mds_bray, province,  label = TRUE, cex=1) # add confidence ellipses around habitat types
  dev.off()
  
#Plot NMDS showing specific factors (region with spider) 
  png(file="NMDS_vegan_Region_spider.png", width=1500, height=1450, res=200)
par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
    mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-1,1), ylim=c(-1, 1))   # don't plot anything yet
    points(mds.fig, "sites", pch = 19, col = "aquamarine", select = Region == 
           "N. Coastal")   # plot just the samples, colour by habitat, pch=19 means plot a circle
    points(mds.fig, "sites", pch = 19, col = "chartreuse", select = Region == 
           "N.gyre")
    points(mds.fig, "sites", pch = 19, col = "royalblue", select = Region == 
           "Equator")   # plot just the samples, colour by habitat, pch=19 means plot a circle
    points(mds.fig, "sites", pch = 19, col = "grey80", select = Region == 
           "S.Coastal")   # plot just the samples, colour by habitat, pch=19 means plot a circle
    points(mds.fig, "sites", pch = 19, col = "tomato", select = Region == 
           "S.gyre")
    ordispider(mds_bray, Region, label = TRUE, cex=1.3) # add confidence ellipses around habitat types
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
  
  ef2<-envfit(mds_bray ~ Florescence, contex, permu=999) # pick environmental factors
  
  png(file="NMDS_Surface_temp_fluor.png", width=1500, height=1450, res=200)
  plot(mds_bray, display="sites",type = "p", xlim=c(-1,1), ylim=c(-1,1))
  points(mds_bray, display = "sites", cex = 0.8, pch=21, col="black", bg="grey80")             
  plot(ef2)
  tmp<-with(contex, ordisurf(mds_bray, Temperature, add=TRUE, col="green")) 
  text(mds_bray, display = "sites", cex=0.7, col="grey20", pos=4)
  dev.off()
  }
  
  results<-list(rank=rank, ef=ef, ef2=ef2)
  return(results)
  
}





NGS_Significance_function= function(x,y,z){

  lat_long <- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1)
    lat_long_dist=vegdist(lat_long, method="euc")
  euk <- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1)
    euk_dist=vegdist( euk, method="jaccard")
  euk_general <- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1)
    euk_general_dist=vegdist(euk_general, method="jaccard")
  euk_alg <- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1)
    euk_alg_dist=vegdist(euk_alg, method="jaccard")
  euk_abun <- read.csv(file.choose(), header = TRUE, sep=";", row.names = 1)
    euk_abun_dist=vegdist(euk_abun, method="jaccard")
  #envrinmental parameter 
  bray_all<- vegdist(comm, method="jaccard")  # calculate Bray-Curtis distance among samples
  mds_bray<<-metaMDS(comm, distance="jaccard") # robustness
ef<= envfit(mds_bray, contex, permu=999) 

anosim(comm_stand, size_fraction, perm=1000, distance="bray")
anosim(comm_stand, province, perm=1000, distance="bray")
anosim(comm_stand, Region, perm=1000, distance="bray")
anosim(comm, Salinity, perm=1000, distance="jaccard")
anosim(comm, Oxy, perm=1000, distance="jaccard")
anosim(comm, Temperature, perm=1000, distance="jaccard")
anosim(comm, C02, perm=1000, distance="jaccard")
anosim(comm, Florescence, perm=1000, distance="jaccard")
anosim(comm, Nitrate, perm=1000, distance="jaccard")
anosim(comm, Phosphate, perm=1000, distance="jaccard")
anosim(comm, Silicate, perm=1000, distance="jaccard")

comm_dist=vegdist(comm_stand, method="bray")
contex_dist=vegdist(contex, method="euc")
mantel(comm_dist, contex_dist, perm=1000) #  Mantel's permutation test for similarity of two matrices
# mantel(comm_dist, euk_general_dist, perm=1000) # also works for two communtiy matrixes

mantel.partial(comm_dist, contex_dist, lat_long_dist, method = "pearson", permutations = 999) # can also be
# dont with accouting for geographic distance.

mantel(comm_dist, euk_general_dist, perm=1000)
mantel(comm_dist, euk_dist, perm=1000)
mantel(comm_dist, euk_alg_dist, perm=1000)
mantel(comm_dist, euk_abun_dist, perm=1000)


}
  
  
  




