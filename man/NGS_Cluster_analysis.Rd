% Generated by roxygen2 (4.1.0.9000): do not edit by hand
% Please edit documentation in R/cluster_analysis.R
\name{NGS_Cluster_analysis}
\alias{NGS_Cluster_analysis}
\title{Cluster analysis}
\usage{
NGS_Cluster_analysis(x)
}
\arguments{
\item{normalised}{community matrix, default na.rm=TRUE}
}
\description{
This function perfroms a dissimilarity analysis using the vegdist function (VEGAN) on a normalised community matrix. It subsequently performs a hierarchical cluster analysis via the (stats) hclust function.
-The function requires a species by site matrix
-It makes labels based on row.names
-Trees are calulated using Jaccard or bray curtis similarity matrixes vegdist()
  and clustered by "single", "average" or "complete" linkage algorithms hclust() functions.
  Trees are then plotted together for comparison.
- The function directly plots the anylsis.
}
\examples{
NGS_Cluster_analysis(x)
}
\keyword{analysis}
\keyword{cluster}

