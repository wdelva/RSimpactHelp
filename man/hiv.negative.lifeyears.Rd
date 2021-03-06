% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hiv.negative.lifeyears.R
\name{hiv.negative.lifeyears}
\alias{hiv.negative.lifeyears}
\title{A function that returns the total number of HIV-negative life years lived
between two points in simulation time, for a particular age group and gender.}
\usage{
hiv.negative.lifeyears(datalist = datalist, agegroup = c(15, 30),
  timewindow = c(15, 30), site = "All")
}
\arguments{
\item{datalist}{The datalist that is produced by \code{\link{readthedata}}}

\item{agegroup}{HIV ngative within this agegroup.}

\item{timewindow}{HIV negative time window e.g timewindow = c(30, 40).}

\item{site}{Select only the particular site from the study, if all ignore site/use all sites.}
}
\value{
the total number of HIV-negative life years lived aggregated by gender.
}
\description{
A function that returns the total number of HIV-negative life years lived
between two points in simulation time, for a particular age group and gender.
}
\examples{
data(datalist)
hiv.negative.lifeyears <- hiv.negative.lifeyears(datalist = datalist,
agegroup=c(15, 40), timewindow=c(15, 30), site="All")
hiv.negative.lifeyears

}
