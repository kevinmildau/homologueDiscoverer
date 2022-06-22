#' sdbCreate
#'
#' Function creates series_db object from annotated_peak_table. Scalar or vector inputs for annotation_id (integer), noise (boolean) and homologue_series_name (character) can be provided to select specific homologue_ids from the annotated peak_table and give them a manual noise status or homologue series name. By default all three inputs are NA, leading to all homologue_ids from the annotated peak table being used, with noise = TRUE and homologue_series_name = NA.
#'
#' @param annotated_peak_table An annotated peak table in tibble format as produced by identifyHomologues.
#' @param annotation_id Integer scalar or vector of ids to include in the seriesdb.
#' @param noise Boolean scalar or vector of same length as annotation id with noise status for homologue_ids in annotation_id.
#' @param homologue_series_name Character scalar or vector of same length as annotation id with names for homologue_ids in annotation_id.
#' @param sample_origin Character string indicating the sample / processing run origin of the identified series in non-standardized.
#' @param sample_description Character string containing any relevant information about the sample / processing run in non-standardized from.
#' @return A series_db tibble.
#' @export
#'
sdbCreate <- function(annotated_peak_table, annotation_id = NA,
                      noise = TRUE, homologue_series_name = NA,
                      sample_origin = NA, sample_description = NA){
  # Check Annotated table information
  checkAnnotatedTableFormat(annotated_peak_table)
  # Remove non-homologue peaks from annotated peak table
  annotated_peak_table <- filter(annotated_peak_table, !(is.na(homologue_id)))
  if (nrow(annotated_peak_table) == 0){
    stop('Error: No non-NA homologue_id values in annotated peak_table.')
  }
  if (is.na(annotation_id)){
    print("No homologue id selection provided. Using all homologue id values from annotated peak table.")
    annotation_id <- annotated_peak_table %>% select(., homologue_id) %>% pull(.) %>% unique(.)
      unique(annotated_peak_table$homologue_id)
  }
  seriesdb <- annotated_peak_table %>%
    # Extract only selected annotation_ids from annotated peak table
    right_join(., tibble(homologue_id = annotation_id), by = "homologue_id") %>%
    # Build and join provided metadata for annotation_ids
    left_join(., tibble(homologue_id = annotation_id,
                        homologue_series_name = homologue_series_name,
                        noise = noise,
                        timestamp = Sys.time(),
                        push_id = as.integer(1),
                        sample_origin = as.character(sample_origin),
                        sample_description = as.character(sample_description)),
              by = "homologue_id") %>%
    mutate(., sample_peak_id = peak_id) %>%
    # Select only relevant seriesdb columns
    select(., mz, rt, intensity, homologue_id, within_series_id, noise,
           homologue_series_name, timestamp, push_id, sample_peak_id, sample_origin,
           sample_description)
  # Find and add max intensity identifier to seriesdb
  max_intensities <- seriesdb %>%
    group_by(., homologue_id) %>% # grouping for max by homologue_id
    slice_max(intensity, n = 1, with_ties = FALSE) %>%
    ungroup() %>% mutate(., max_intensity = TRUE) %>%
    select(., sample_peak_id, homologue_id, max_intensity)
  # Join max identifier with seriesdb & add normalized intensities
  seriesdb <- full_join(seriesdb, max_intensities,
                        by = c("homologue_id", "sample_peak_id")) %>%
    mutate(., max_intensity = if_else(is.na(max_intensity), FALSE, max_intensity)) %>%
    group_by(., homologue_id) %>%
    mutate(., normalized_intensity = intensity / sum(intensity)) %>%
    select(., -"intensity") %>%
    ungroup(.) %>%
    arrange(., homologue_id, within_series_id)
  if(nrow(seriesdb) == 0){
    stop("Input conditions lead to empty seriesdb. Database creation aborted.")
  }
  return(seriesdb)
}
