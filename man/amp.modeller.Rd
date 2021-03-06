% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amp.modeller.R
\name{amp.modeller}
\alias{amp.modeller}
\title{Model age-mixing pattern.}
\usage{
amp.modeller(dataframe, agegroup, timepoint, timewindow, start = FALSE,
  SHIMS = TRUE, method = "lmer")
}
\arguments{
\item{dataframe}{The dataframe that is produced by
\code{\link{agemix.episodes.df.maker}}}

\item{agegroup}{Boundaries of the age group that should be retained, e.g.
c(15, 30). The interval is closed on the left and open on the right.}

\item{timepoint}{Point in time during the simulation to be used in the
calculation.}

\item{timewindow}{The length of time before the timepoint for which
relationships should be included, e.g. 1, representing one year before the
timepoint. This should be a whole number.}

\item{start}{This is a logical indicating that only relationships starting
after the beginning of the window should be used. If start = FALSE
relationships could start before the time window. This is the default.}

\item{SHIMS}{Should only the three most recent relationships be kept, like in
the SHIMS survey?}

\item{method}{Should the \code{link[nlme]{lme}} or \code{link[lme4]{lmer}}
function be used? The lme function can be used to model heteroskedastic
residual error.}
}
\value{
returns a model fitted with \code{link[nlme]{lme}} or
  \code{link[lme4]{lmer}}
}
\description{
Models the age-mixing pattern for the population simulated in Simpact. The
user specifies a time point, time window, and age group for which they would
like to obtain a summary of the age-mixing pattern. The user also specifies
whether only relationships that started within the time window should be
used.
}
\details{
The \code{amp.modeller} function can use either the \code{link[nlme]{lme}} or
the \code{link[lme4]{lmer}} function to build a linear mixed effects model
regressing the partner's age at the time the relationship started on his/her
own age at the beginning of the relationship. The models are stratified by
gender. Each person can have more than one relationship so a random intercept
at the level of the person is added to the model.

The models built with \code{link[nlme]{lme}} also explicitly account for
heteroskedastic variance, since the variance in partner ages may grow for
older ages.

The function produces a fitted model.
}
\examples{
data(persreldf)
agemixpatdat <- amp.modeller(dataframe = persreldf, agegroup = c(18, 50),
timepoint = 30, timewindow = 0.5, start = FALSE, SHIMS = TRUE, method = "lmer")

}
