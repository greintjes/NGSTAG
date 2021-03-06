% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/species_occurance.R
\name{NGS_species_occurance_function}
\alias{NGS_species_occurance_function}
\title{Species_occurance_function}
\usage{
NGS_species_occurance_function(x, echo = T)
}
\arguments{
\item{normalised}{community matrix, echo=TRUE}
}
\description{
This function anaylsis a normalised community matrix.
- It takes a species by site matrix and analyses species occurance within and between sites.
- Optional transformation of the data matrix to emphasise or deemphasis dominants
- How many sites does each species occur in ? This part of the function first calls (veg=0)
  which evaluates True/flase or 1/0, and then sums the column values
- Species occurance all sites or 1 or less sites (spc.pres X, X being number of sites)
- Gets the average cover for each species (sum in site divided by the number of sites)
- optional identify which are the dominant organism using listpts function
- Calculate species area relationship
- all results are stored in file speceis_occurance_results.csv includes names of species in
  all sites or in 1 or less sites. Also top 10 abundance species are included. All graphs
  (.png images) are stored in the working directory. They are high quailty images
}
\examples{
NGS_species_occurance_functiom(x, echo=TRUE)
}
\keyword{occurance}
\keyword{species}

