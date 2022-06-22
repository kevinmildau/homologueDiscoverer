#' sdbSummarize
#'
#' Function generates a summary table with information about the different homologues in seriesdb. The summary table is printed with 10 decimals as well as returned as a tibble.
#'
#' For each homologue id in seriesdb, the summary contains the series_length (number of peaks in series), the mean and median differences in mz between peaks, as well as the molecular formula estimate and the theoretical increment estimate as computed by MassTools::CalcMF.
#'
#' @param seriesdb A seriesdb tibble created with createSeriesDB().
#' @param calcMF_tolerance Maximum ppm tolerance for MassTools::CalcMF.
#'
#' @return Summary table in tibble format.
#' @export
sdbSummarize <- function(seriesdb, calcMF_tolerance = 100){
  old_print_settings <- options(pillar.sigfig = 10)
  summary <- seriesdb %>%
    group_by(., homologue_id) %>%
    summarize(., series_length = n(),
              mean_diff = mean(diff(mz)),
              median_diff = median(diff(mz)),
              mass_tools =
                wrapperCalcMF(mean_diff,
                                      tolerance = calcMF_tolerance)) %>%
    unpack(., mass_tools)
  print(summary, n = Inf)
  options(old_print_settings)
  return(summary)
}

#' wrapperCalcMF
#'
#' Helper function for sdbSummarize that uses MassTools::calcMF on a provided mz value with given tolerance, and extracts the hit with the lowest ppm, if any. Results are returned as a tibble with one row and a molecular formula as well as a theoretical increment.
#'
#' @param mz Mass to charge value for which to calculate Molecular Formula and Theoretical Mass.
#' @param tolerance Maximum ppm tolerance for MassTools::CalcMF.
#'
#' @return Tibble with one row and two columns for molecular_formula (character) and theoretical_increment (double).
#' @keywords internal
wrapperCalcMF <- function(mz, tolerance){
  df <- MassTools::calcMF(mz, z = 0, tolerance)
  if (!(is.null(df))){
    df <- df %>%
      mutate(., ppm = abs(ppm)) %>%
      slice_min(., ppm)
    return(tibble(molecular_formula = select(df, MF) %>% pull(),
                  theoretical_increment = select(df, mz) %>% pull()))
  } else {
    return(tibble(molecular_formula = as.character(NA),
                  theoretical_increment = as.double(NA)))
  }
}
