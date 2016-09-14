#' Fraction of women within specified age group with more than one sexual partner in the past year.
#'
#' \code{degree.df.frac()} Computes the fraction of women within specified age
#' group with more than one partner from the degreedata dataframe.
#'
#' @param degreedata.df The dataframe that is produced by \code{\link{degree.df.maker()}}
#' @param survey.time Time point of the cross-sectional survey.
#' @param agegroup Boundaries of the age group (lower bound <= age < upper bound) that should be retained, e.g. c(15, 30)
#' @param window.width Time period before the survey e.g 1 year before the survey.
#' @param only.new Logical indicator. If TRUE, only relationships that were newly started during window.width are counted
#' (i.e. the individual was NEVER in a relationship with these partners before the start of the window).
#' If FALSE, all relationships that were ongoing at some point during the window.width are counted.
#'
#' @return A value, that is, the fraction of women with >1 sexual partner in the last year.
#'
#' @examples
#' load(degreedata.df)
#' degreefrac.df <- degree.df.frac(degreedata.df, agegroup = c(15, 30), survey.time = 10, window.width = 1, only.new = TRUE)

degree.df.frac<- function(degreedata.df){
  # ntotal<-dim(degree.df)[1]
  ntotal<-dim(degreedata.df)[1]  # all women aged 15-30 in relationships in the last year.
  degree2plus.df<-filter(degreedata.df, Degree>=2)
  npartial<-dim(degree2plus.df)[1]    # women aged 15-30 with more than 1 sexual partner in the last year.
  degree.frac<-npartial/ntotal   # fraction of women aged 15-30 with >1 partner in the last year.
  return(degree.frac)
}

