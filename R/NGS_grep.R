#' NGSgrep
#'
#' This function allows you to grep a specific string from a column name. For example a specific taxnomic name or groups. All the columns containing the string and their associated values are copied into a new vector. 
#' 
#' @param x, vector
#' @keywords grep
#' @export 
#' @examples
#' NGSgrep(x)
NGSgrep<-function(x){
greep=grep("string", names(x), value = TRUE)

  greeped=x[, c(greep)]
  greeped=data.matrix(greeped)
write.table(greeped, "string.csv", sep=";")
}