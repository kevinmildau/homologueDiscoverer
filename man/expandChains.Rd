% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_detectHomologues.R
\name{expandChains}
\alias{expandChains}
\title{expandChains}
\usage{
expandChains(
  chain_list,
  peak_table,
  tolerance,
  rel_rt_tolerance_high,
  rel_rt_tolerance_low
)
}
\arguments{
\item{chain_list}{A list of data frames / tibbles containing the the chains to be extended provided appropriate candidate peaks are present in peak table.}

\item{peak_table}{A matrix version of the peak table with mz (double), rt (double), and peak_id (int) columns.}

\item{tolerance}{Error tolerance in ppm (error relative to last chain node mz + mz stepsize).}

\item{rel_rt_tolerance_high}{Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates, the value is used as a symmetric margin. If the chain is length three, rel_rt_tolerance_high is used in conjunction with rel_rt_tolerance_low to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing; static rt step sizes are accounted for using rel_rt_tolerance_low). Used in findChainCandidates.}

\item{rel_rt_tolerance_low}{Value between 0 and 1. Sets the rt tolerance relative to the current rt step of a explored chain. If the chain has only two candidates it is ignored. If the chain is length three, rel_rt_tolerance_low is used in conjunction with rel_rt_tolerance_high to create asymmetric retention time margins based on the initialized series rt step trends (monotonically increasing or decreasing). rel_rt_tolerance_low sets the rt margin in the opposite direction of the series trend, thus allowing for static rt steps or mild deviation from monotonic trends. Used in findChainCandidates.}
}
\value{
list of chain data frames / tibbles. In general, this will contain at least one data frame (initialized chain of two elements). Can contain a larger number of data frames if many branching events possible.
}
\description{
Helper function that runs the main recursive element the package on initialized chains.
}
\keyword{internal}
