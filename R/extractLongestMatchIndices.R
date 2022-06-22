#' extractLongestMatchIndices
#'
#' Helper Functions evaluates a series of mz values with respect to their repetitive increment. Returns the longest sequence of mz values that follow each other within tolerance of the step size. If two sequences of mz values matching the stepsize criterion are of equal length, the first one is returned.
#'
#' @param mz_vals Vector of mz values.
#' @param min_match_length Minimum match length.
#' @param ppm_tolerance Tolerance in ppm for step size between any two mz vals.
#' @param step_size Mz step size to evaluate diff(mz_vals) against.
#'
#' @return Indices of mz_vals for which meet the longest subsequence criterion. If none, a boolean FALSE is returned.
#' @keywords internal
extractLongestMatchIndices <- function(mz_vals, min_match_length, ppm_tolerance, step_size){
  if (length(mz_vals) == 0 | length(mz_vals) < min_match_length){
    return(FALSE)
  }
  differences <- diff(mz_vals)
  mz_ref <- max(mz_vals, na.rm = TRUE) # ppm_tolerance to mz tolerance uses most liberal reference mz
  mz_tolerance <- getMzTolerance(mz_ref, ppm_tolerance)
  boolout <- near(step_size, differences, tol = mz_tolerance)
  rleboolout <- rle(boolout)
  lengths <- rleboolout$lengths
  values <- rleboolout$values
  max_seq_bool <- with(rleboolout,
                       rep(lengths == max(lengths[values]) &
                             values &
                             lengths >= min_match_length, lengths))
  # There is a possibility for two max length subsequences being found;
  if (!any(max_seq_bool == FALSE)){
    return(seq(1:length(mz_vals))) # all diffs within tolerance, no gaps
  } else if (!any(max_seq_bool)){
    # no differences in expected range, no sequence found OR
    # no sequence of differences equal long enough for a match
    return(FALSE)
  } else {
    first_true <- min(which(max_seq_bool == TRUE)) # find first element of first max-seq
    subset_bool <- max_seq_bool[seq(first_true,length(max_seq_bool))]
    if (!any(subset_bool == FALSE)){
      # Only initial mismatches
      return (seq(first_true, length(mz_vals)))
    } else {
      switch_point <- first_true + min(which(max_seq_bool[first_true:length(max_seq_bool)] == FALSE))
      first_longest_seq_indices <- c(first_true, switch_point - 1)
      mz_indices <- seq(first_longest_seq_indices[1], first_longest_seq_indices[2])
      return(mz_indices)
    }
  }
}
