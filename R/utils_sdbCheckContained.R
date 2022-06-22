#' checkAnyMatch
#'
#' Helper function checks whether there is any match within ppm and rt tolerance of some target mz and rt tuple within the target df. Target df needs to contain two columns: mz, rt. Function returns NA if no match, and the homologue_id corresponding to the mz-rt tuple if there is one or more matches.
#'
#' @param target_mz Target mass to charge ratio.
#' @param target_rt Target retention time.
#' @param homologue_id Current peak homologue id.
#' @param ppm_tolerance Tolerance in part per million.
#' @param rt_tolerance Retention time tolerance. Beware of peak table retention time units!
#' @param peak_table A peak table (tibble) with mz and rt columns at least.
#'
#' @return NA or the matched homologue_id
#' @keywords internal
checkAnyMatch <- function(target_mz, target_rt, homologue_id, ppm_tolerance, rt_tolerance, peak_table){
  mz_tolerance <- getMzTolerance(target_mz, ppm_tolerance)
  match <- any(near(peak_table[, "rt"], target_rt, tol = rt_tolerance) &
                 near(peak_table[, "mz"], target_mz, tol = mz_tolerance))
  if (match) {
    return(homologue_id)
  } else {
    return(as.integer(NA))
  }
}
