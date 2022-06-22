#' annotateHomologues
#'
#' Function loops through homologue_ids in seriesd_db and checks whether they are contained within the peak table. If contained, the peak table rows are annotated with the corresponding homologue_id from series_db.
#'
#' annotateHomologues performs overlap checks in multiple stages. First, for any peaks in series_db corresponding to a homologue_id of interest, the closest match in peak_table abiding by ppm_tolerance and rt_tolerance is extracted.
#' Once all peaks were matched, a second evaluation process extracts the longest gap-free alignment within the matching peaks abiding by a ppm_tolerance in the step size separating each peak. Annotation is only done if there is a gap-free alignment larger than min_match_length. In order to prevent annotation corruption through overwriting, annotation is only allowed if no previous annotation exists for the matched peaks. Multiple series annotations for the same peak may occur as a result of too relaxed matching criteria or through series_db objects composed of multiple runs for which no overlap checks were done.
#'
#' A quick_search toggle has been introduced to indicate whether a maximum intensity peak pre-screening of the peak table should be performed. The idea behind this heuristic is that if a homologue within series_db is to be contained within the peak table, then it's maximum intensity peak from the series_db should also be present. If not present, the alignment step is skipped for the homologue_id in question.
#'
#' @param peak_table A peak table with mz, rt, intensity and peak_id columns. Intensity is not used by the algorithm but is needed for downstream functionality. Additional columns are allowed.
#' @param seriesd_db A series_db tibble created with sdbCreate().
#' @param ppm_matching_tolerance Tolerance in ppm for matching two peaks.
#' @param ppm_step_tolerance Tolerance in ppm relative to largest mz value in homologue series for evaluating matched peaks increment consistency.
#' @param rt_tolerance Retention time tolerance in matching. Beware of rt units in the database!
#' @param min_match_length Minimum gap free overlap between two homologue series for them to be considered overlapping. Should be >= 3.
#' @param quick_search Boolean toggle indicating whether search should be made faster by pre-screening the max intensity peaks of the seriesd_db against the peak_table. If the highest intensity peak of the homologue in series_db is not present in the peak table, then the series itself is also unlikely to be present, allowing for a reduction of the number of peak matching tasks.
#'
#' @return Peak table as a tibble with peaks annotated with homologue_identifiers from series_db using a homologue_id and additional homologue information. If no overlap detected, the input peak_table is returned without modifications.
#' @export
annotateHomologues <- function(peak_table, seriesd_db, ppm_matching_tolerance,
                             ppm_step_tolerance, rt_tolerance,
                             min_match_length, quick_search = TRUE){
  # Augment Peak table to contain row_id & annotation column
  peak_table <- peak_table %>%
    rowid_to_column(., var = "ptb_row_id") %>% # add row_id for annotation purposes
    mutate(annotation = as.integer(NA),
           sdb_row_id = as.integer(NA))
  # Augment series_db to contain sdb_row_id & sdb_sample_peak_id columns.
  seriesd_db <- seriesd_db %>%
    rowid_to_column(., var = "sdb_row_id") %>%
    mutate(., sdb_sample_peak_id = sample_peak_id)
  # Pre-screen peak table if quick search is TRUE
  if (quick_search){
    seriesd_db <- preScreenPeakTable(peak_table, seriesd_db, ppm_matching_tolerance, rt_tolerance)
  }
  # If any series left after pre screening, perform indepth screening.
  if (nrow(seriesd_db) > 0){
    series_ids <- unique(seriesd_db$homologue_id)
    for (series in series_ids){
      # Extract homologue series data from series_dbdb
      series_data <- seriesd_db[seriesd_db[,"homologue_id"] == series,, drop = FALSE]
      # find all matches in peak_table and return their row_id and mz values
      matches <- mapply(findClosestMatch, series_data$mz, series_data$rt, series_data$sdb_row_id,
                        MoreArgs = list(ppm_matching_tolerance, rt_tolerance, peak_table),
                        SIMPLIFY = FALSE) %>%
        bind_rows(.) %>% na.omit(.)
      diff <- mean(series_data %>% pull(mz) %>% diff(.)) # mean difference of the series homologue.
      # Extract longest continuous (gap free) series of matching diff in the homologue series, if any
      indices <- extractLongestMatchIndices(mz_vals = pull(matches, mz),
                                            min_match_length, ppm_step_tolerance,
                                            diff)
      if (!is.logical(indices)){
        matches <- matches[indices, ] # extract the longest matching subsequence of matches
        ptb_indices <- pull(matches, "ptb_row_id")
        # peak_table[ptb_indices, c("annotation")] %>% pull()
        annotation_tmp <- peak_table[ptb_indices, c("annotation")] %>% pull()
        if (!(all(is.na(annotation_tmp)))){
          warning(paste("Attempted annotation of already annotated peaks in peak table. Annotation of series_id", series, "in peak table is prevented. Overlap concerns peak_table rows:", toString(ptb_indices), ". Already Present annotations are the following:", toString(annotation_tmp), ". This issue can be caused by single peaks having multiple annotations in the series database, or peak matching settings being too relaxed, allowing for spurious matching."))
        } else {
          peak_table[pull(matches, "ptb_row_id"), c("annotation", "sdb_row_id")] <-
            cbind(series, pull(matches, sdb_row_id))
        }
      }
    }
    # Check whether any annotations were made
    if (all(is.na(peak_table$annotation))){
      print("No overlap found in screening. No peaks annotated to belong to a homologue. Returning input peak table.")
      return(select(peak_table, -c("ptb_row_id", "annotation")))
    } else {
      annotated_peak_table <-
        as_tibble(peak_table) %>%
        arrange(., annotation, mz) %>%
        mutate(., homologue_id = annotation) %>%
        mutate(., homologue_id = if_else(is.na(homologue_id), 0L, homologue_id)) %>%
        left_join(., select(seriesd_db, sdb_row_id, within_series_id, noise,
                            homologue_series_name, sdb_sample_peak_id), by = c("sdb_row_id")) %>%
        select(., -c("sdb_row_id", "ptb_row_id")) %>%
        mutate(., noise = if_else(is.na(noise), FALSE, noise))
      return(annotated_peak_table)
    }
  } else {
    print("No overlap found in pre-screening. No peaks annotated to belong to a homologue. Returning input peak table.")
    return(select(peak_table, -c("sdb_row_id", "annotation")))
  }
}

