% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/NGS_merge.R
\name{NGSmerge}
\alias{NGSmerge}
\title{Merging function}
\usage{
NGSmerge(x, y)
}
\arguments{
\item{c_merge_taxonomy(x,y,)}{X =species site vector 1, y= species site vector 2. The specific function takes two vectors and combines their taxonomy and make a new vector called merge1. this can be done multiple times to add several vector togehter.}

\item{d_sort.removeNA(x),}{X= merge vector with NA. NA are replaced with 0's and the vector is ordered.}
}
\description{
This function allows you to merge several dataframes based on their taxonomy (colname). This function is based on the rbind.all.columns() from http://amywhiteheadresearch.wordpress.com/2013/05/13/combining-dataframes-when-the-columns-dont-match/. It finds the difference between the column names between two vectors (setdiff()) and then adds dummy columns to each vector representing the missing columns and adds NA into the column as data. Then it uses rbind() to combines the vectors.
}
\keyword{NGSmerge(x,y)}
\keyword{merge}

