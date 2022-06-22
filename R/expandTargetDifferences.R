#' expandTargetDifferences
#'
#' Allows for a provided set of mass differences to be expanded to account for possible higher charges of the molecules.
#'
#' When looking for an increment of 14, an increment of 7 may also be relevant if the molecule has charge 2. This function creates a list of all possible increments provided theoretical mass increments and a vector of charges to be considered. 14 and 1,2,3 -> 14, 7, 4.667
#'
#' @param diff_mz A list of theoretical mass increments.
#' @param charges A list of charges ranging from 1 to n, where n is the maximum charge considered.
#'
#' @return List of unique increments to be used in the targeted homologue search.
#' @export
#'
#' @examples
#' expandTargetDifferences(14, c(1,2,3))
#'
expandTargetDifferences <- function(diff_mz, charges){
  increments <- expand_grid(diff_mz, charges) %>%
    mutate(., increment = diff_mz / charges) %>%
    distinct(., increment) %>%
    pull(increment)
  return(increments)
}
