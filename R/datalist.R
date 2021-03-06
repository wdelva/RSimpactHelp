#' Simpact run datalist data.
#'
#' This is an example datalist that could be created from simpact.run.
#' The interventions are enabled as documented in the standard simulation
#' as defined in \code{\link{intervention.introduced}} and the age distribution
#'  are standard as set in \code{\link{agedistr.creator}}.
#' The simulation is run on 1000 men and 1000 women for 40 years on seed 2.
#'
#' @docType data
#'
#' @usage data(datalist)
#'
#'
#' @format  The datalist has 7 elements, each are being described below.
#' \describe{
#' \item{ptable}{Person Table with 21 columns and  2242 rows. This keeps
#' information about each individual in the simulatoin with the simulation window.}
#' \item{rtable}{Relationship table with 6 columns and 535 rows. Each relationship
#' formed and desolved is registered in this table. With the age difference being
#' calculated between the two partners.}
#' \item{etable}{All the events within the simulation that took place are registered
#' in this table. 10 columns and 4488 rows.}
#' \item{ttable}{Treatment table}
#' \item{itable}{Information about the simulation configuration parameters}
#' \item{ltable}{Population log and the number of people on treatment}
#' \item{vltable}{Viral Load log for each person whose viral load is subjected
#' to change and the description on why that change is happening}
#' }

"datalist"

