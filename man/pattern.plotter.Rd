% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pattern.plotter.R
\name{pattern.plotter}
\alias{pattern.plotter}
\title{Plot age-mixing pattern.}
\usage{
pattern.plotter(dl)
}
\arguments{
\item{dl}{The datalist that is produced by \code{\link{pattern.modeller}}}
}
\value{
a large ggplot object of the age-mixing pattern
}
\description{
Visualizes the age-mixing pattern from a Simpact-simulated population. This
function takes as input a datalist produced by \code{\link{pattern.modeller}}
and outputs a scatterplot of the age of the individual in the population
versus the age of their partner. The plot contains a dashed line that
represents a scenario where people choose partners that are the same age as
themselves. Juxtaposed with this is a line that represents the
population-average predicted partner ages based upon the generalised linear
mixed effects model in \code{\link{pattern.modeller}}.
}
