
#' Subset sequences based on gender ratio and age of individuals
#'
#' @param simpact.trans.net Transmission networks computed by \code{\link{transmission.network.builder()}}
#' @param limitTransmEvents Choose transmission network with at least this  amount of individuals
#' @param timewindow Time interval in which sampling are considered
#' @param seq.cov Sequece coverage
#' @param age.limit Age limit for all individuals
#' @return A vector of selected IDs
#' @importFrom dplyr filter
#' @export
#'


IDs.Seq.Random <- function(simpact.trans.net = simpact.trans.net,
                           limitTransmEvents = 7,
                           timewindow = c(10,40),
                           seq.cov = 70,
                           age.limit=65){

  seeds.id <- length(simpact.trans.net)

  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)

  for(i in 1:seeds.id){

    transm.age.i <- as.data.frame(simpact.trans.net[[i]])

    age.i <- transm.age.i$SampTime - transm.age.i$TOBRec

    transm.age.i <- cbind(transm.age.i, age.i)

    new.transm.tab[[i]] <- transm.age.i

  }

  # ID numbers of Selected networks with at least limitTransmEvents + 1 indiviuals

  IDs.transm <- vector()

  TransmEventsCountVector <- vector()

  for(k in 1:seeds.id){
    trans.net.i.check <- as.data.frame(new.transm.tab[[k]])

    if(nrow(trans.net.i.check)>=limitTransmEvents){

      TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))

      IDs.transm <- c(IDs.transm, k)
    }
  }

  if(length(IDs.transm)>=1){

    ## Binding together all selected transmission transmission networks ##

    for (q in 1:length(IDs.transm)){

      if(q==1){
        p <- IDs.transm[q]
        trans.sum <- new.transm.tab[[p]]
        rename.id <- paste0(p,".",trans.sum$id,".C")
        trans.sum$id <- rename.id
        trans.sum.rename.id <- trans.sum
      }
      else{

        p <- IDs.transm[q]

        read.trans.sum <- new.transm.tab[[p]]
        rename.id.read <- paste0(p,".",read.trans.sum$id,".C")
        read.trans.sum$id <- rename.id.read
        trans.sum.rename.id <- rbind(trans.sum.rename.id, read.trans.sum)
      }

    }

    trans.sum.age.limit <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$age.i<=age.limit)

    trans.sum.men <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$SampTime >= timewindow[1] & trans.sum.age.limit$SampTime <= timewindow[2])

    trans.sum.women <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$SampTime >= timewindow[1] & trans.sum.age.limit$SampTime <= timewindow[2])

    perc.100 <- nrow(trans.sum.men) + nrow(trans.sum.women) # total number of individuals with age limit

    trans.sum.men.women <- rbind(trans.sum.men, trans.sum.women)

    perc.seq.coverage <- round(perc.100*seq.cov/100) # total number of wanted individuals at seq.cov sequence coverage


    samp.all <- sample(trans.sum.men.women$id, perc.seq.coverage)

  }else{
    samp.all <- NA
  }

  return(samp.all)
}

# IDs.Seq.Random merge

