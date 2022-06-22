#' detectHomologues
#'
#' Main function that searches a provided peak_table for homologues in the specified mode. A peak table (tibble) is returned with each repetitive series annotated by an integer identifier, NA representing no series membership (all non-homologue).
#'
#'  Main function for homologue detection. Defines outermost while loop over peak table.
#'
#'  Algorithm outline:
#'
#'  While there are peaks in the ordered peak table
#'
#'  --> Use the first indexed peak and provided search parameters to initialize peak chain 2-tuples.
#'
#'  --> If valid chain initialization:
#'
#'  ---> Do recursive chain expansion. Using search parameters, expand 2-tuples recursively until no further candidates match the search restrictions.
#'
#'  ---> Of all chains created, select longest
#'
#'  ---> Eval chain length above min length, if so, declare homologue
#'
#'  ---> Remove first indexed peak or all members of the homologue series and repeat.
#'
#' Homologue series heuristics are built into candidate search, hence no further checks for homologue series trend consistency are done.
#' @param peak_table A peak table (type tibble) with mz (double), rt (double), and peak_id (int) columns.
#' @param mz_min Minimum mz increment for a new peak to be considered (untargeted, absolute value). Must be 0 or larger. Needs to be set for untargeted search mode.
#' @param mz_max Maximum mz increment for a new peak to be considered (untargeted, absolute value). Must be mz_min or larger. Needs to be set for untargeted search mode.
#' @param rt_min Minimum rt increment for a new peak to be considered. Must be 0 or larger.
#' @param rt_max Maximum rt increment for a new peak to be considered. Must be rt_min or larger.
#' @param ppm_tolerance Error tolerance in ppm (error relative to last chain node mz + mz step-size).
#' @param search_mode "targeted" or "untargeted" determining the search approach.
#' @param mz_steps Vector of mz step sizes to be searched around. Needs to be set for targeted search mode.
#' @param min_series_length Minimum size of a repetitive series for it to be considered a homologue. The minimum plausible value is 3. Any value below 5 is likely to return spurious results unless ppm constrains are highly rigid.
#' @param step_mode "increment" or "decrement" to indicate whether mz steps are to be looked for in an increasing or decreasing mass-to-charge ratio value over retention time fashion.
#' @param rel_rt_tolerance_high Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates, the value is used as a symmetric margin. If the chain is length three, rel_rt_tolerance_high is used in conjunction with rel_rt_tolerance_low to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing; static rt step sizes are accounted for using rel_rt_tolerance_low). Used in findChainCandidates.
#' @param rel_rt_tolerance_low Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates it is ignored. If the chain is length three, rel_rt_tolerance_low is used in conjunction with rel_rt_tolerance_high to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing). rel_rt_tolerance_low sets the rt margin in the opposite direction of the series trend, thus allowing for static rt steps or mild deviation from monotonic trends. Used in findChainCandidates.
#' @param progress_counter Controls console update frequency of remaining peaks messages. Default is 100.
#' @param verbose Can be used to deactivate progress counter messages in console. Default is TRUE.
#' @return Peak table with additional homologue_id column with integer identifiers for each identified series. An homologue_id of NA indicating no series membership.
#' @export
#'
#' @examples
#' # Targeted Run:
#' detectHomologues(peak_table, mz_min = NA, mz_max = NA, rt_min, rt_max, ppm_tolerance, search_mode = "targeted", mz_steps = c(14.0156, 44.02628), min_series_length = 5, step_mode = "increment")
#' # Untargeted Run:
#' detectHomologues(peak_table, mz_min = 10, mz_max = 20, rt_min, rt_max, ppm_tolerance, search_mode = "untargeted", mz_steps = NA, min_series_length = 5, step_mode = "increment")
detectHomologues <- function(peak_table,
                             mz_min = NA,
                             mz_max = NA,
                             rt_min,
                             rt_max,
                             ppm_tolerance,
                             search_mode = NA,
                             mz_steps = NA,
                             min_series_length = 4,
                             step_mode = "increment",
                             rel_rt_tolerance_high = 0.6,
                             rel_rt_tolerance_low = 0.1,
                             progress_counter = 100,
                             verbose = T){
  # Input checks
  checkPeakTableFormat(peak_table)
  checkDetectHomologueSearchParameters(mz_min, mz_max, rt_min,
                                       rt_max, ppm_tolerance,
                                       search_mode, mz_steps,
                                       min_series_length, step_mode)
  # Initializing search
  if (step_mode == "increment") {
    peak_table <- peak_table %>% arrange(., mz,rt)
  } else if (step_mode == "decrement"){
    peak_table <- peak_table %>% arrange(., desc(mz),rt)
  }
  annotated <- peak_table %>% mutate(homologue_id = as.integer(NA))
  peak_table <- matrix(c(peak_table$mz, peak_table$rt, peak_table$peak_id), ncol = 3)
  colnames(peak_table) <- c("mz", "rt", "peak_id")
  n_homologues = 0
  while (dim(peak_table)[1] >= 1){
    chain_list <- initializeChains(peak_table, mz_min, mz_max, rt_min, rt_max,
                                   ppm_tolerance, search_mode, mz_steps,
                                   step_mode)
    if (length(chain_list) != 0){
      chain_list <- expandChains(chain_list, peak_table,
                                 ppm_tolerance, rel_rt_tolerance_high,
                                 rel_rt_tolerance_low)
      selected_chain <- extractLongestChain(chain_list)
      if (nrow(selected_chain) >= min_series_length){
        if(verbose){print("Found Homologue Series!")}
        tmp_ids <- selected_chain[, "peak_id"]
        peak_table <-
          peak_table[!(peak_table[,"peak_id"] %in% tmp_ids),,drop = FALSE]
        n_homologues = as.integer(n_homologues + 1)
        annotated <- annotated %>%
          mutate(homologue_id = if_else(peak_id %in% tmp_ids, n_homologues,
                                      homologue_id))
        if(verbose){printProgress(peak_table, progress_counter)}
      } else {
        peak_table <- peak_table[-1,,drop=FALSE]
        if(verbose){printProgress(peak_table, progress_counter)}
      }
    } else {
      peak_table <- peak_table[-1,,drop=FALSE]
      if(verbose){printProgress(peak_table, progress_counter)}
    }
  }
  annotated <- annotated %>%
    group_by(., homologue_id) %>%
    mutate(., within_series_id = row_number()) %>%
    mutate(., within_series_id = if_else(homologue_id == 0L, as.integer(NA),
                                         within_series_id)) %>%
    ungroup(.)
  return(annotated)
}


