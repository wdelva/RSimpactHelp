#' Read the .csv files and combine them into one list
#'
#' @param modeloutput The object that is produced by \code{\link{RSimpactCyan::simpact.run()}}
#' @return A list, containing dataframes for the output of the model run:
#' ptable (people), rtable (relationships), etable (events), ttable (HIV treatment episodes),
#' paramtable (input parameters), and ltable (bookkeeping log) if created
#' @examples
#' cfg <- initiate()
#' modeloutput <- simpact.run(configParams = cfg, destDir = "~TMP")
#' datalist <- readthedata(modeloutput)


readthedata <- function(modeloutput){
  path <- as.character(modeloutput["outputfile"])
  outputID <- as.character(modeloutput["id"])
  DestDir <- sub(pattern = paste0(outputID, "output.txt"), replacement = "", x = path, fixed=T)
  personlogfilename <- paste0(DestDir, outputID, "personlog.csv")
  relationlogfilename <- paste0(DestDir, outputID, "relationlog.csv")
  eventlogfilename <- paste0(DestDir, outputID, "eventlog.csv")
  treatmentlogfilename <- paste0(DestDir, outputID, "treatmentlog.csv")
  inputparamlogfilename <- paste0(DestDir, outputID, "settingslog.csv")
  periodiclogfilename <- paste0(DestDir, outputID, "periodiclog.csv")

  ptable <- data.table::fread(personlogfilename, sep = ",", skip = 0)
  rtable <- data.table::fread(relationlogfilename, sep = ",", skip = 0)
  readetable <- readcsvcolumns::read.csv.columns(eventlogfilename, has.header = FALSE, column.types = "rssiirsiir")
  etable <- data.table::setDT(readetable)
  etable.colnames <- colnames(etable)
  if (ncol(etable) == 10){
    data.table::setnames(etable, etable.colnames,
             c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))
  } else {
    data.table::setnames(etable, etable.colnames,
             c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age", "extradescript", "value"))

  }
  ttable <- data.table::fread(treatmentlogfilename, sep  = ",", skip = 0)
  itable <- data.table::fread(inputparamlogfilename, sep  = ",", skip = 0)

  if (file.exists(periodiclogfilename)){
    ltable <- data.table::fread(periodiclogfilename, sep = ",", skip = 0)
    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, itable = itable, ltable = ltable)
  } else {
    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, itable = itable)
  }
  return(outputtables)

}


degree.df.maker <- function(datalist = datalist, survey.time = 10, window.width = 1, only.new = TRUE){
  # For the people in the datalist,
  # make a dataframe that has the variables of the ptable,
  # but in addition, also the number of relationships that were ongoing (if only.new == FALSE)
  # or that were newly formed (if only.new == TRUE) at some time in the survey window.
  # (E.g. The one-year period from 9 to 10 years into the simulation).

  return(degree.df)
}

degree.hist.maker <- function(dataframe = degree.df){
  # Some plotting commands, e.g. using ggplot2,
  # to visualise the degree distribution.
}

degree.distr.fitter <- function(dataframe = degree.df){
  # Fit various parametric distribution to the degree data
  # and evaluate / compare their goodness-of-fit.
  # Normal, Poisson, Negative Bionomial, Gamma, Beta, Log-normal, Weibull, ...
  # And return the parameters (e.g. mean and variance or scale and shape) of the
  # best fitting distribution.
}




