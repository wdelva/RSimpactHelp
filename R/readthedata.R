#' Read the Simpact output files.
#'
#' Read the .csv files and combine them into one list.
#'
#' @param modeloutput The object that is produced by RSimpactCyan::simpact.run.
#' @return A list, containing dataframes for the output of the model run:
#' ptable (people), rtable (relationships), etable (events), ttable (HIV treatment episodes),
#' itable (input parameters), ltable (bookkeeping log) if created and vltable (HIV viral load)
#' @examples
#' cfg <- list()
#' modeloutput <- RSimpactCyan::simpact.run(configParams = cfg, destDir = "temp")
#' datalist <- readthedata(modeloutput)
#'
#' @importFrom data.table fread
#' @importFrom data.table setDT
#' @importFrom data.table setnames
#' @importFrom readcsvcolumns read.csv.columns
#' @export

readthedata <- function(modeloutput){
  
  path <- as.character(modeloutput["outputfile"])
  outputID <- as.character(modeloutput["id"])
  DestDir <- sub(pattern = paste0(outputID, "output.txt"), 
                 replacement = "", 
                 x = path, 
                 fixed = T)
  personfile <- paste0(DestDir, outputID, "personlog.csv")
  relationfile <- paste0(DestDir, outputID, "relationlog.csv")
  eventfile <- paste0(DestDir, outputID, "eventlog.csv")
  treatmentfile <- paste0(DestDir, outputID, "treatmentlog.csv")
  inputparamfile <- paste0(DestDir, outputID, "settingslog.csv")
  periodicfile <- paste0(DestDir, outputID, "periodiclog.csv")
  viralloadfile <- paste0(DestDir, outputID, "hivviralloadlog.csv")
  facilitiesxyfile <- paste0(DestDir, outputID, "facilitypositions.csv")

  ptable <- data.table::fread(personfile, sep = ",", skip = 0)
  vltable <- data.table::fread(viralloadfile, sep = ",", skip = 0)
  rtable <- data.table::fread(relationfile, sep = ",", skip = 0)
  readetable <- readcsvcolumns::read.csv.columns(eventfile, 
                                                 has.header = FALSE,
                                                 column.types = "rssiirsiir")
  etable <- data.table::setDT(readetable)
  etable.colnames <- colnames(etable)
  
  if (ncol(etable) == 10){
    data.table::setnames(etable, etable.colnames,
             c("eventtime", "eventname", "p1name", "p1ID", "p1gender",
               "p1age", "p2name", "p2ID", "p2gender", "p2age"))
  } else {
    data.table::setnames(etable, etable.colnames,
             c("eventtime", "eventname", "p1name", "p1ID", "p1gender",
               "p1age", "p2name", "p2ID", "p2gender", "p2age", 
               "extradescript", "value"))

  }
  
  ttable <- data.table::fread(treatmentfile, sep  = ",", skip = 0)
  itable <- data.table::fread(inputparamfile, sep  = ",", skip = 0)

  if (file.exists(periodicfile) && file.exists(facilitiesxyfile)){
    
    ltable <- data.table::fread(periodicfile, sep = ",", skip = 0)

    ftable <- data.table::fread(facilitiesxyfile, sep = ",", skip = 0)

    outputtables <- list(ptable = ptable, 
                         rtable = rtable, 
                         etable = etable,
                         ttable = ttable, 
                         itable = itable, 
                         ltable = ltable,
                         vltable = vltable, 
                         ftable = ftable)
    
  }else if(file.exists(facilitiesxyfile) && !file.exists(periodicfile)){
    
    ftable <- data.table::fread(facilitiesxyfile, sep = ",", skip = 0)

    outputtables <- list(ptable = ptable, 
                         rtable = rtable, 
                         etable = etable,
                         ttable = ttable, 
                         itable = itable, 
                         vltable = vltable, 
                         ftable = ftable)

  }else if(!file.exists(facilitiesxyfile) && file.exists(periodicfile)){
    
    ltable <- data.table::fread(periodicfile, sep = ",", skip = 0)

    outputtables <- list(ptable = ptable, 
                         rtable = rtable, 
                         etable = etable,
                         ttable = ttable, 
                         itable = itable, 
                         ltable = ltable, 
                         vltable = vltable)
  }else{
    outputtables <- list(ptable = ptable, 
                         rtable = rtable,
                         etable = etable, 
                         ttable = ttable, 
                         itable = itable, 
                         vltable = vltable)
  }

  return(outputtables)
}
