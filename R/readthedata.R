#' Read the Simpact output files.
#'
#' Read the .csv files and combine them into one list.
#'
#' @param modeloutput The object that is produced by \code{\link{RSimpactCyan::simpact.run()}}
#' @return A list, containing dataframes for the output of the model run:
#' ptable (people), rtable (relationships), etable (events), ttable (HIV treatment episodes),
#' itable (input parameters), ltable (bookkeeping log) if created and vltable (HIV viral load)
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
  viralloadlogfilename <- paste0(DestDir, outputID, "hivviralloadlog.csv")

  ptable <- data.table::fread(personlogfilename, sep = ",", skip = 0)
  vltable <- data.table::fread(viralloadlogfilename, sep = ",", skip = 0)
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

  data.table::setnames(vltable, vltable.colnames, c("VLTimeLog", "ID", "Description","Log10SPVL", "Log10VL"))

  if (file.exists(periodiclogfilename)){
    ltable <- data.table::fread(periodiclogfilename, sep = ",", skip = 0)
    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, itable = itable, ltable = ltable, vltable = vltable)
  } else {
    outputtables <- list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable, itable = itable, vltable = vltable)
  }
  return(outputtables)

}
