% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database_management_functions.R
\name{sdbSave}
\alias{sdbSave}
\title{sdbSave}
\usage{
sdbSave(series_db, filename, check_exists = TRUE)
}
\arguments{
\item{series_db}{A series_db tibble as generated by sdbCreate()}

\item{filename}{A absolute or relative file path specifying the location and name of the new series_db file. Filename should end on .RDS}

\item{check_exists}{Toggle for file existence check defaults to TRUE. Set to FALSE for automatic file replacement without warning.}
}
\value{
None. series_db object is written to hard drive.
}
\description{
Function saves series_db to file using provided absolute or relative file path as the name. Function performs existence checks to warn user of possible overwriting of existing files using interactive prompt.
}
\details{
The interactive prompt can be deactivated using check_exists. See function input parameters.
}
