#' Normalisation of community data
#'
#' This function can be used to normalise both community and contextual data using the vegan decostand () function. Choose the method of normailsation wisely !!
#' Its exports the normalised community matrix and a normalisation "checked" matrix.
#' @param input .csv community or environmental matrix
#' @keywords normalisation
#' @export
#' @examples
#' NGS_normalise_abundance_function(x)
#' 

NGS_normalise_abundance_function =function(x){  
  community_stand<<- decostand(x, method = "total")  # "total" divide by sample total abundance, "norm" make margin sum of 
  # squares equal to one
  norm_check<<-apply(community_stand, 1, sum)            # check total abundance in each sample, should be all 1 now.
  
  
  
}
