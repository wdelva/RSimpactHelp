% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ART.coverage.calculator.R
\name{ART.coverage.calculator}
\alias{ART.coverage.calculator}
\title{Calculate overall HIV prevalence and ART coverage aggregated by gender.}
\usage{
ART.coverage.calculator(datalist = datalist, agegroup = c(15, 30),
  timepoint = 30, site = "All")
}
\arguments{
\item{datalist}{The datalist that is produced by \code{\link{readthedata}}}

\item{agegroup}{Boundaries of the age group (lower.bound <= age < upper.bound) that
should be retained, e.g. agegroup = c(15, 30)}

\item{timepoint}{Point in time at which the ART coverage should be calculated.}

\item{site}{Select only the particular site from the study, if all ignore site/use all sites.}
}
\value{
a dataframe with HIV prevalence estimates and ART coverage estimate and surrounding
confidence bounds,
for the specified time point and age group, overall, and stratified by gender
}
\description{
Calculate the HIV prevalence and ART coverage at a point in time, for specific
age groups and gender strata.
}
\examples{
data(datalist)
ART.coverage.df <- ART.coverage.calculator(datalist = datalist,
agegroup = c(15, 30), timepoint = 30, site="All")

}
