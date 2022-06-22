#' findClosestMatch
#'
#' Helper function that finds closest match within rt and ppm tolerance of some target mz and rt tuple in a peak_table. Extracts the row_id and mz value for further processing.
#'
#' @param target_mz Target mz value around which to search.
#' @param target_rt Target rt value around which to search.
#' @param sdb_row_id Seriesdb row identifier.
#' @param ppm_tolerance MZ tolerance for matching in ppm.
#' @param rt_tolerance Retention time tolerance for matching.
#' @param peak_table Peak table within which to search.
#'
#' @return A tibble with 1 row and three columns; peak table row id, sdb row id and mz value.
#' @keywords internal
findClosestMatch <- function(target_mz, target_rt, sdb_row_id, ppm_tolerance, rt_tolerance, peak_table){
  mz_tolerance <- getMzTolerance(target_mz, ppm_tolerance)
  matches_within_tolerance <-
    peak_table[ near(peak_table[, "rt"], target_rt, tol = rt_tolerance) &
                  near(peak_table[, "mz"], target_mz, tol = mz_tolerance),,
                drop = FALSE ]
  if (nrow(matches_within_tolerance) != 0) {
    # Extract closest mz match within tolerance paramaters
    abs_diff = abs( (select(matches_within_tolerance, mz) %>% pull())  - target_mz)
    match <- matches_within_tolerance[which.min(abs_diff),, drop = FALSE]
    match <- select(match, ptb_row_id, mz) %>%
      rename(., ptb_row_id = ptb_row_id) %>%
      mutate(., sdb_row_id = sdb_row_id)
    return(match)
  } else {
    return(tibble("ptb_row_id" = as.integer(NA),
                  "mz" = as.double(NA),
                  "sdb_row_id" = as.integer(NA)))
  }
}



