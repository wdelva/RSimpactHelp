#' Fraction of women within specified age group with more than one
#' sexual partner in the past year.
#'
#' \code{degree.df.frac} Computes the fraction of women within specified age
#' group with more than one partner from the degreedata dataframe.
#'
#' @param degreedata.df The dataframe that is produced by \code{\link{degree.df.maker}}
#' @return A value, that is, the fraction of women with >1 sexual partner in the last year.
#'
#@examples
#data(degreedata.df)
#degreefrac.df <- degree.df.frac(degreedata.df)
#' @import dplyr
#' @export

degree.df.frac<- function(degreedata.df){
  # all women aged 15-30 in relationships in the last year. ???
  ntotal <- dim(degreedata.df)[1]
  degree2plus.df <- dplyr::filter(degreedata.df, Degree>=2)

  #women aged 15-30 with more than 1sexual partner in the last yr.
  npartial <- dim(degree2plus.df)[1]

  # fraction of women aged 15-30 with >1 partner in the last yr.
  degree.frac <- npartial/ntotal

  return(degree.frac)
}

# This example is not correct
#degreefrac.df <- degree.df.frac(degreedata.df, agegroup = c(15, 30), survey.time = 10,
#window.width = 1, only.new = TRUE)
