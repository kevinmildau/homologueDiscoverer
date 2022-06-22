#' getMzTolerance
#'
#' Helper function that computes the absolute mz tolerance given a mz and tolerance in ppm.
#'
#' @param mz MZ value for which the tolerance margin is to be computed.
#' @param ppm_tolerance Tolerance in ppm.
#'
#' @return Absolute tolerance in mz.
#' @keywords internal
getMzTolerance <- function(mz, ppm_tolerance){
  mz_tolerance <- ((mz * (1 + (ppm_tolerance / 10^6))) -
                     (mz * (1 - (ppm_tolerance / 10^6))) ) / 2
  return(mz_tolerance)
}
