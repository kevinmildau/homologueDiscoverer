#' isScalar
#'
#' Function checks whether input is atomic, and if so, of length 1.
#'
#' @param x Any R object.
#'
#' @return Boolean True or False.
#' @keywords internal
isScalar <- function(x){
  return(is.atomic(x) && length(x) == 1L)
}
