% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client.facility.R
\name{client.facility}
\alias{client.facility}
\title{A function that determines the facility the client attended
Assumption is that the client facility will never change once asssigned}
\usage{
client.facility(datalist = datalist, site = "All")
}
\arguments{
\item{datalist}{The datalist that is produced by \code{\link{readthedata}}}
}
\value{
the datalist with facility attached to each client
}
\description{
TO REVISE (Need to get the first time that ART was started then reset when
drop out and look at the VL six month after, Ignore everthing else
after this)
}
\examples{
data(datalist)
datalist <- client.facility(datalist = datalist, site="All")
datalist

}
