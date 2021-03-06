% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ART.coverage.vector.creator.R
\name{ART.coverage.vector.creator}
\alias{ART.coverage.vector.creator}
\title{Create a vector of ART coverage measurements.}
\usage{
ART.coverage.vector.creator(datalist = datalist, agegroup = c(15, 50))
}
\arguments{
\item{datalist}{The datalist that is produced by \code{\link{readthedata}}}

\item{agegroup}{Boundaries of the age group (lower bound <= age < upper
bound) that should be retained, e.g. c(15, 30)}
}
\value{
A vector with ART coverage estimates, for the specific age group.
}
\description{
Calculate the ART coverage, for a specific age group, in yearly intervals,
from the start of the simulation until the end of the simulation.
}
\examples{
data(datalist)
art.cov.vector.df <- ART.coverage.vector.creator(datalist = datalist, agegroup = c(15,50))
art.cov.vector.df

}
