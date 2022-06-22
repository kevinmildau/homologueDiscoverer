#' checkAnnotatedTableFormat
#'
#' @param annotated_table peak table with annotation column as produced by detectHomologues.
#'
#' @return TRUE if annotated table matches assertions. Stops program if not.
#' @keywords internal
checkAnnotatedTableFormat <- function(annotated_table){
  tibble_condition <- is_tibble(annotated_table)
  nrow_condition <- nrow(annotated_table) > 0
  column_condition <- all(c("mz", "rt", "homologue_id", "within_series_id",
                            "intensity", "peak_id") %in%
                            names(annotated_table))
  stopifnot("Provided table must be of type tibble." = tibble_condition)
  stopifnot("Provided table must have nrow > 0." = nrow_condition)
  stopifnot('All of the following column names must be in provided table: "mz", "rt", "homologue_id", "within_series_id", "intensity", "peak_id"' = column_condition )
  return(TRUE)
}



#' checkSeriesDBFormat
#'
#' Helper function checks whether seriesdb input corresponds to format produced by createSeriesDB.
#'
#' @param seriesdb seriesdb object as produced by sdbCreate()
#'
#' @return TRUE if seriesdb matches format assumptions. Stops program if not.
#' @keywords internal
checkSeriesDBFormat <- function(seriesdb){
  tibble_condition <- is_tibble(seriesdb)
  nrow_condition <- nrow(seriesdb) > 0
  column_condition <- all(c("mz", "rt", "homologue_id",
                            "within_series_id", "noise", "homologue_series_name",
                            "timestamp", "push_id", "max_intensity", "normalized_intensity",
                            "sample_peak_id", "sample_origin",
                            "sample_description") %in%
                            names(seriesdb))
  stopifnot("Provided table must be of type tibble." = tibble_condition)
  stopifnot("Provided table must have nrow > 0." = nrow_condition)
  stopifnot('All of the following column names must be in provided table: "mz", "rt", "homologue_id", "within_series_id", "noise", "homologue_series_name", "timestamp", "push_id", "max_intensity", "normalized_intensity", "sample_peak_id", "sample_origin", "sample_description"' = column_condition )
  return(TRUE)
}
