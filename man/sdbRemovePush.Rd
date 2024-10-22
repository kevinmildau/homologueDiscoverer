% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database_management_functions.R
\name{sdbRemovePush}
\alias{sdbRemovePush}
\title{sdbRemovePush}
\usage{
sdbRemovePush(series_db, push_ids)
}
\arguments{
\item{series_db}{A series_db tibble as generated by sdbCreate()}

\item{push_ids}{A single integer or list/vector of push_id integers to be filtered from series_db.}
}
\value{
Filtered series_db tibble. If filtering removes all rows of series_db, NA will be returned alongside a warning message.
}
\description{
Function filters series_db object to not contain any of the push_ids provided (scalar or list/vector of integers).
}
