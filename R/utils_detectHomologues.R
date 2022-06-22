#' nearSelect
#'
#' Helper function filters down a peak table (matrix format) to contain only mz rows falling within upper bounds and lower bounds determined via the target mz value and an error tolerance in ppm.
#'
#' Context: is used in vectorized nearSelect. Used in initializeChains.
#'
#' @param value The target mz value for searching.
#' @param table Peak table in matrix form, with named column "mz".
#' @param tolerance Relative error tolerance specified in ppm.
#'
#' @return Table filtered down to those rows matching the search criteria. May be empty, nrow = 0, but remains a data frame type object.
#' @keywords internal
#'
nearSelect <- function(value, table, tolerance){
  # value is the target mz value. tolerance is the acceptable error in ppm
  LB <- value * (1 - (tolerance / 10^6)) # lower bound using ppm tolerance
  UB <- value * (1 + (tolerance / 10^6)) # upper bound using ppm tolerance
  return(table[table[,"mz"] >= LB & table[,"mz"] <= UB,,drop = FALSE])
}



#' printProgress
#'
#' Helper function which prints progress dependent on the number of rows left in peak table (matrix) and a progress counter indicting the message frequency.
#'
#' @param peak_table Peak table in matrix form.
#' @param progress_counter Integer value determining progress print. If nrow peak table is a multiple of progress counter, the remaining number of peaks is printed.
#'
#' @keywords internal
printProgress <- function(peak_table, progress_counter){
  if((nrow(peak_table) %% progress_counter) == 0){print(paste("Remaining Peaks:", nrow(peak_table)))}
}



#' customBindRows
#'
#' Helper function that allows vectorized uses of bindRows over lists of chain candidates to be combined together. Each chain candidate row i from table y is merged with an existing table x. Candidates are rows of y, the index i is given by a lapply call. The result of repeated uses of this function is a list of tables x, with each x being rowbound to a row of y.
#'
#' Context: Function is used in initializeChains.
#'
#' @param i Row index as provided by a lapply call of a index list.
#' @param x Data frame / tibble of peaks (existing chain).
#' @param y Data frame / tibble of peaks (chain candidates to be added)
#'
#' @return Row bound version of x and a row of y.
#' @keywords internal
#'
customBindRows <- function(i, x, y){
  return(rbind(x, y[i,,drop=FALSE]))
}



#' initializeTargeted
#'
#' Function creates a (likely small) number of candidate series tables to be explored in expandChains().
#'
#' The targeted approach creates a very narrow window using the target increments and ppm tolerances, as well as a rt tolerance to find candidate peaks to form seconds steps with respect to the root peak. Creates an initial chain for each pairwise combination of root peak + any node contained within these narrow windows.
#'
#' @param peak_table A matrix version of the peak table with mz (double), rt (double), and peak_id (int) columns.
#' @param wrt_min Minimum rt increment for a new peak to be considered.
#' @param wrt_max Maximum rt increment for a new peak to be considered.
#' @param tolerance Error tolerance in ppm (error relative to last chain node mz + mz stepsize).
#' @param target_increments Vector of increments to be searched around (targeted)
#' @param mode2 "increment" or "decrement" to indicate whether mz steps are to be looked for in an incremental or decremental way.
#'
#' @return list of data frames containing each two rows that are a possible series start. Empty list if no candidate matches.
#' @keywords internal
#'
initializeTargeted <- function(peak_table, wrt_min, wrt_max, tolerance, target_increments, mode2){
  root <- peak_table[1,, drop=FALSE] # take new root node
  if (mode2 == "increment"){
    candidate_mz <- root[, "mz"] + target_increments # define the target mz values to look for (with error tolerance)
  } else if (mode2 == "decrement"){
    candidate_mz <- root[, "mz"] - target_increments # define the target mz values to look for (with error tolerance)
  } else {
    stop("No valid mode2 setting provided. Please provide 'increment' or 'decrement'.")
  }
  rtfiltered_peak_table <- peak_table[peak_table[,"rt"] >= root[,"rt"] + wrt_min &
                                        peak_table[,"rt"] <= root[,"rt"] + wrt_max, ,drop=FALSE] # rt filtering
  stepll <- lapply(candidate_mz, nearSelect, rtfiltered_peak_table, tolerance) # near select with tolerance
  stepll <- base::Filter(function(x){nrow(x) >= 1}, stepll) # this removes any empty matrices
  # this is a multirange check, where each range check may introduce another matrix of hits,
  if (length(stepll) > 1){
    steptb <- do.call(rbind, stepll)
    chain_list <- lapply(1:nrow(steptb), customBindRows, root, steptb)
    return(chain_list)
  } else  if (length(stepll) == 1) {
    if (nrow(stepll[[1]] > 1)){
      steptb <- do.call(rbind, stepll)
      chain_list <- lapply(1:nrow(steptb), customBindRows, root, steptb)
    } else {
      chain_list <- list(rbind(root, stepll[[1]]))
      return(chain_list)
    }
  } else {
    return(list())
  }
}



#' initializeUntargeted
#'
#' Function creates a (possibly large) number of candidate series tables to be explored in expandChains(). The approach creates a filter window using mz and rt bounds, and creates an initial chain for each pairwise combination of root + any node in these windows.
#'
#' @param peak_table A matrix version of the peak table with mz (double), rt (double), and peak_id (int) columns.
#' @param wmz_min Minimum mz increment for a new peak to be considered (untargeted, absolute value)
#' @param wmz_max Maximum mz increment for a new peak to be considered (absolute value, untargeted)
#' @param wrt_min Minimum rt increment for a new peak to be considered.
#' @param wrt_max Maximum rt increment for a new peak to be considered.
#' @param tolerance Error tolerance in ppm (error relative to last chain node mz + mz stepsize).
#' @param mode2 "increment" or "decrement" to indicate whether mz steps are to be looked for in an incremental or decremental way.
#'
#' @return list of data frames containing each two rows that are a possible series start. Empty list if no candidate matches.
#' @keywords internal
initializeUntargeted <- function(peak_table, wmz_min, wmz_max, wrt_min, wrt_max, tolerance, mode2){
  root <- peak_table[1,, drop=FALSE] # take new root node
  if (mode2 == "increment"){
    filtered_peak_table <- peak_table[peak_table[,"rt"] >= root[,"rt"] + wrt_min &
                                        peak_table[,"rt"] <= root[,"rt"] + wrt_max &
                                        peak_table[,"mz"] >= root[,"mz"] + wmz_min &
                                        peak_table[,"mz"] <= root[,"mz"] + wmz_max, ,drop=FALSE]
  } else if (mode2 == "decrement"){
    filtered_peak_table <- peak_table[peak_table[,"rt"] >= root[,"rt"] + wrt_min &
                                        peak_table[,"rt"] <= root[,"rt"] + wrt_max &
                                        peak_table[,"mz"] <= root[,"mz"] - wmz_min &
                                        peak_table[,"mz"] >= root[,"mz"] - wmz_max, ,drop=FALSE]
  } else {
    stop("No valid mode2 setting provided. Please provide 'increment' or 'decrement'.")
  }
  if (nrow(filtered_peak_table) > 1){
    chain_list <- lapply(1:nrow(filtered_peak_table), customBindRows, root, filtered_peak_table)
    return(chain_list)
  } else  if (nrow(filtered_peak_table) == 1) {
    chain_list <- list(rbind(root, filtered_peak_table))
    return(chain_list)
  } else {
    return(list())
  }
}



#' initializeChains
#'
#' Function creates initial chain tables for the current root node (first row of peak_table). Two modes are available: "targeted" or "untargeted". In both cases, each root may have more than one possible candidate second element. Hence, a list of possible chains is returned to be followed in subsequent recursive calls of expandChains.
#' @param peak_table A matrix version of the peak table with mz (double), rt (double), and peak_id (int) columns.
#' @param wmz_min Minimum mz increment for a new peak to be considered (untargeted, absolute value)
#' @param wmz_max Maximum mz increment for a new peak to be considered (absolute value, untargeted)
#' @param wrt_min Minimum rt increment for a new peak to be considered.
#' @param wrt_max Maximum rt increment for a new peak to be considered.
#' @param tolerance Error tolerance in ppm (error relative to last chain node mz + mz stepsize).
#' @param mode "targeted" or "untargeted" determining the search approach.
#' @param target_increments Vector of increments to be searched around (targeted)
#' @param mode2 "increment" or "decrement" to indicate whether mz steps are to be looked for in an incremental or decremental way.
#'
#' @return list of data frames containing each two rows that are a possible series start. Empty list if no candidate matches.
#' @keywords internal
initializeChains <- function(peak_table, wmz_min, wmz_max, wrt_min, wrt_max,
                             tolerance, mode, target_increments, mode2){
  # mode can be : "targeted" or "untargeted"
  if (!(mode %in% c("targeted", "untargeted"))){
    stop("Please select a valid mode, options are: 'targeted' and 'untargeted'.")
  }
  if (mode == "targeted"){
    return(initializeTargeted(peak_table, wrt_min, wrt_max,
                              tolerance, target_increments, mode2))
  }
  if (mode == "untargeted"){
    return(initializeUntargeted(peak_table, wmz_min, wmz_max, wrt_min, wrt_max,
                                tolerance, mode2))
  }
}



#' expandChains
#'
#' Helper function that runs the main recursive element the package on initialized chains.
#'
#' @param chain_list A list of data frames / tibbles containing the the chains to be extended provided appropriate candidate peaks are present in peak table.
#' @param peak_table A matrix version of the peak table with mz (double), rt (double), and peak_id (int) columns.
#' @param tolerance Error tolerance in ppm (error relative to last chain node mz + mz stepsize).
#' @param rel_rt_tolerance_high Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates, the value is used as a symmetric margin. If the chain is length three, rel_rt_tolerance_high is used in conjunction with rel_rt_tolerance_low to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing; static rt step sizes are accounted for using rel_rt_tolerance_low). Used in findChainCandidates.
#' @param rel_rt_tolerance_low Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates it is ignored. If the chain is length three, rel_rt_tolerance_low is used in conjunction with rel_rt_tolerance_high to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing). rel_rt_tolerance_low sets the rt margin in the opposite direction of the series trend, thus allowing for static rt steps or mild deviation from monotonic trends. Used in findChainCandidates.
#'
#' @return list of chain data frames / tibbles. In general, this will contain at least one data frame (initialized chain of two elements). Can contain a larger number of data frames if many branching events possible.
#'
#' @keywords internal
expandChains <- function(chain_list, peak_table, tolerance, rel_rt_tolerance_high, rel_rt_tolerance_low){
  # Recursive function that grows a list of chains with leave nodes a long as possible.
  # Only the longest possible chain will be kept (leave) by means of the outer most merge call;
  # as long as a chain can be grown, it will be followed and grown.
  tmp <- list()
  for (chaintb in chain_list){
    # finding chains for an initialized chain works differently from initialization of a chain, hence another function
    # a tibble is returned here, if empty, there are no chain candidates,
    chain_candidates <- findChainCandidates(chaintb, peak_table, tolerance, rel_rt_tolerance_high, rel_rt_tolerance_low)
    if (length(chain_candidates) != 0) {
      tmp <- c(tmp, expandChains(chain_candidates, peak_table, tolerance, rel_rt_tolerance_high, rel_rt_tolerance_low)) # recursive call
    } else {
      tmp <- append(tmp, list(chaintb)) # the end node tibble is kept.
    }
  }
  return(tmp)
}



#' findChainCandidates
#'
#' Helper function that extends an chain data frame of two or more peaks (root, first step, second step, etc.) provided there are suitable candidate peaks in the peak table.
#'
#' Context: Helper function within expandChain that checks for new candidates in peak table at each chain extension step. Reqursion starts after chain initialization with expandChain and this findChainCandidates.
#'
#' Empty list return is processed in expandChains as no chain extension, keeping the originally supplied data frame. This is done since returning the input data frame would require additional length checks to validate that the chain was grown.
#'
#' @param chaintb An initialized chain data frame.
#' @param peak_table A matrix version of the peak table with mz (double), rt (double), and peak_id (int) columns.
#' @param tolerance Error tolerance in ppm (error relative to last chain node mz + mz stepsize).
#' @param rel_rt_tolerance_high Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates, the value is used as a symmetric margin. If the chain is length three, rel_rt_tolerance_high is used in conjunction with rel_rt_tolerance_low to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing; static rt step sizes are accounted for using rel_rt_tolerance_low). Used in findChainCandidates.
#' @param rel_rt_tolerance_low Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates it is ignored. If the chain is length three, rel_rt_tolerance_low is used in conjunction with rel_rt_tolerance_high to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing). rel_rt_tolerance_low sets the rt margin in the opposite direction of the series trend, thus allowing for static rt steps or mild deviation from monotonic trends. Used in findChainCandidates.
#' @return List of tibbles containing extended chain data frames. If no candidates, returns empty list.
#' @keywords internal
findChainCandidates <- function(chaintb, peak_table, tolerance,
                                rel_rt_tolerance_high, rel_rt_tolerance_low){
  # Function within recursive findLongestChain; part of lvl 3 to n.
  # Finds continuation candidates for the current peak chain and returns
  # them as a tibble.

  # Get last peak in chain mz data
  current_mz <- chaintb[nrow(chaintb), "mz"]
  current_rt <- chaintb[nrow(chaintb), "rt"]
  # Get difference in rt between last two entries (based on latest entries in chain)
  rt_step <- abs(chaintb[nrow(chaintb)-1, "rt"] - chaintb[nrow(chaintb), "rt"])
  mz_step <- chaintb[, "mz"]%>% diff() %>% mean(.) # mean observed step size among accepted series members

  if (nrow(chaintb) >= 3) {
    # Enforce Chain rt continuity (monotonically increasing or monotonically decreasing)
    # Based on first encountered trend between step 1 and step 2
    rt_step_1 <- abs(chaintb[1, "rt"] - chaintb[2, "rt"])
    rt_step_2 <- abs(chaintb[2, "rt"] - chaintb[3, "rt"])
    rt_step_trend = rt_step_1 - rt_step_2 # positive if decreasing, negative if increasing

    if (rt_step_trend >= 0) {
      # only look for rt_step + decreasing margin relative to rtstep
      chain_candidates <- peak_table[peak_table[, "mz"] >= (current_mz + mz_step) * (1 - (tolerance / 10^6)) &
                                       peak_table[, "mz"] <= (current_mz + mz_step) * (1 + (tolerance / 10^6)) &
                                       peak_table[, "rt"] >= current_rt + rt_step - rel_rt_tolerance_high * rt_step &
                                       peak_table[, "rt"] <= current_rt + rt_step + rel_rt_tolerance_low * rt_step,, drop = FALSE]
    } else {
      # only look for rt_step + increasing margin relative to rt_step;
      chain_candidates <- peak_table[peak_table[, "mz"] >= (current_mz + mz_step) * (1 - (tolerance / 10^6)) &
                                       peak_table[, "mz"] <= (current_mz + mz_step) * (1 + (tolerance / 10^6)) &
                                       peak_table[, "rt"] >= current_rt + rt_step - rel_rt_tolerance_low * rt_step &
                                       peak_table[, "rt"] <= current_rt + rt_step + rel_rt_tolerance_high * rt_step,, drop = FALSE]
    }
  } else {
    # Use more lenient rt_step criterion allowing for large positive and negative margin around rt after first step
    chain_candidates <- peak_table[peak_table[, "mz"] >= (current_mz + mz_step) * (1 - (tolerance / 10^6)) &
                                     peak_table[, "mz"] <= (current_mz + mz_step) * (1 + (tolerance / 10^6)) &
                                     peak_table[, "rt"] >= current_rt + rt_step - rel_rt_tolerance_high * rt_step &
                                     peak_table[, "rt"] <= current_rt + rt_step + rel_rt_tolerance_high * rt_step,, drop = FALSE]
  }
  # Output Pre-processing
  n_candidates <- nrow(chain_candidates)
  if(n_candidates > 1){
    return(lapply(1:n_candidates, customBindRows, chaintb, chain_candidates))
  } else if (n_candidates == 1){
    return(list(rbind(chaintb, chain_candidates[1,,drop=FALSE])))
  }else {
    return (list())
  }
}



#' extractLongestChain
#'
#' Helper function that extracts the data frame (data frame or tibble) with the largest number of rows from a list of data frames.
#' In case of a draw the first indexed chain with maximum number of rows is extracted.
#'
#' Context: This function helps resolve multiple chain recursions by selecting the longest running series of peaks.
#'
#' @param chain_list A list of data frames or tibbles, each representing a potential homologue series.
#'
#' @return A data frame or tibble representing a potential homologue series.
#'@keywords internal
extractLongestChain <- function(chain_list){
  chaintb <- chain_list[which.max(sapply(chain_list, nrow))][[1]]
  return(chaintb)
}



#' checkPeakTableFormat
#'
#' Helper function that checks whether the provided peak table
#'
#' @param peak_table
#'
#' @return Boolean True if peak table has the right format.
#' @keywords internal
checkPeakTableFormat <- function(peak_table){
  tibble_condition <- is_tibble(peak_table)
  nrow_condition <- nrow(peak_table) > 0
  column_condition <- all(c("mz", "rt", "intensity", "peak_id") %in%
                            names(peak_table))
  stopifnot("Provided table must be of type tibble." = tibble_condition)
  stopifnot("Provided table must have nrow > 0." = nrow_condition)
  stopifnot('All of the following column names must be in provided table: "mz", "rt", "intensity", "peak_id"' = column_condition )
  return(TRUE)
}



#' checkDetectHomologueSearchParameters
#'
#' Helper function that does rough checks for input parameter validity.
#'
#' @param mz_min Minimum mz increment for a new peak to be considered (untargeted, absolute value). Must be 0 or larger. Needs to be set for untargeted search mode.
#' @param mz_max Maximum mz increment for a new peak to be considered (untargeted, absolute value). Must be mz_min or larger. Needs to be set for untargeted search mode.
#' @param rt_min Minimum rt increment for a new peak to be considered. Must be 0 or larger.
#' @param rt_max Maximum rt increment for a new peak to be considered. Must be rt_min or larger.
#' @param ppm_tolerance Error tolerance in ppm (error relative to last chain node mz + mz step-size).
#' @param search_mode "targeted" or "untargeted" determining the search approach.
#' @param mz_steps Vector of mz step sizes to be searched around. Needs to be set for targeted search mode.
#' @param min_series_length Minimum size of a repetitive series for it to be considered a homologue. The minimum plausible value is 3. Any value below 5 is likely to return spurious results unless ppm constrains are highly rigid.
#' @param step_mode "increment" or "decrement" to indicate whether mz steps are to be looked for in an increasing or decreasing mass-to-charge ratio value over retention time fashion.
#'
#' @return Boolean TRUE if all settings in valid ranges.
#' @keywords internal
checkDetectHomologueSearchParameters <- function(mz_min, mz_max, rt_min,
                                                 rt_max, ppm_tolerance,
                                                 search_mode, mz_steps,
                                                 min_series_length, step_mode){
  stopifnot("rt_min must be scalar!" = isScalar(rt_min))
  stopifnot("rt_max must be scalar!" = isScalar(rt_max))
  stopifnot("mz_min must be scalar!" = isScalar(mz_min))
  stopifnot("mz_max must be scalar!" = isScalar(mz_max))
  stopifnot("ppm_tolerance must be scalar!" = isScalar(ppm_tolerance))
  stopifnot("min_series_length must be scalar!" = isScalar(min_series_length))
  stopifnot("rt_max must be larger than rt_min" = rt_max > rt_min)
  stopifnot("rt_min must be larger than 0" = rt_min > 0)
  stopifnot("min_series_length must be larger than 3!" = min_series_length > 3)
  stopifnot("NA search_mode or step_mode detected. Modes must be defined!" = !is.na(search_mode) & !is.na(step_mode))
  stopifnot("ppm_tolerance must be larger than 0!" = ppm_tolerance > 0)
  if (search_mode == "targeted"){
    stopifnot("mz_steps must be provided for targeted search mode!" = !is.logical(mz_steps))
    if (!is.logical(mz_min) | !is.logical(mz_min)){
      warning("mz ranges provided, yet targted search_mode selected. Ranges will be ignored.")
    }
  } else if (search_mode == "untargeted"){
    stopifnot("mz_min must be larger than 0" = mz_min > 0)
    stopifnot("mz_max must be larger than mz_min" = mz_max > mz_min)
    stopifnot("mz range parameters mz_min and mz_max must be provided for untargeted search mode!" = !is.na(mz_min) & !is.na(mz_max))
    if (!is.logical(mz_steps)){
      warning("mz_steps provided, yet untargeted search_mode selected. Step sizes will be ignored.")
    }
  }
}
