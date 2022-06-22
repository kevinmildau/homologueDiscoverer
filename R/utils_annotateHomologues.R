#' preScreenPeakTable
#'
#' Helper function of annotateHomologues that checks whether seriesdb max intensity peaks of the various homologues have a match within ppm and rt tolerance in the peak table. Filters seriesdb down to those homologue ids which do have a match for complete series matching in annotateHomologues.
#'
#'
#' @param peak_table A peak table with mz and rt columnds.
#' @param seriesdb A seriesdb tibble created with createSeriesDB().
#' @param ppm_tolerance Tolerance in mz for matching two peaks.
#' @param rt_tolerance Retention time tolerance in matching. Beware of rt units in the database!
#'
#' @keywords internal
preScreenPeakTable <- function(peak_table, seriesdb, ppm_tolerance, rt_tolerance){
  # Function checks whether seriesdb homologue id median, max intensity peak is in
  # peak table. If so, the table id is searched more thoroughly, if not, it is
  # skipped from the thorough search to save time.
  # Returns: filtered seriesdb

  # Extract max intensity peaks from seriesdb
  seriesdb_max <- seriesdb %>%
    filter(., max_intensity == TRUE)

  # Lapply a near select for each peak. Return true or false based on whether a match is found.
  matched_ids <- seriesdb_max %>%
    rowwise() %>%
    mutate(matchFound =
             any(near(peak_table[, "rt"], rt, tol = rt_tolerance) &
                   near(peak_table[, "mz"], mz,
                        tol = getMzTolerance(mz, ppm_tolerance)))) %>%
    filter(., matchFound == TRUE) %>%
    select(., homologue_id) %>% pull()
  seriesdb <- seriesdb %>% filter(., homologue_id %in% matched_ids)
  return(seriesdb)
}
