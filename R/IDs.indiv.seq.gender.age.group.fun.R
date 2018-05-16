
#' Subset sequences based on gender ratio and age groups of individuals
#'
#' @param simpact.trans.net Transmission networks computed by \code{\link{transmission.network.builder()}}
#' @param limitTransmEvents Choose transmission network with at least this  amount of individuals
#' @param perc.men Precentage of men to consider when selecting the sequences
#' @param seq.cov Sequece coverage
#' @param age.men Age group for men
#' @param age.women Age group for women
#' @return A list with a vector of ID's of selected sequences, number of men, and women for the selcted sequences and the ratio between men and women among the selcted sequences because in some settings we may have less or more female or male individuals
#' @importFrom dplyr filter
#' @export
#'


IDs.indiv.seq.gender.age.group.fun <- function(simpact.trans.net = simpact.trans.net,
                                               limitTransmEvents = 3,
                                               perc.men = 30,
                                               seq.cov = 35,
                                               age.men=c(15,60),
                                               age.women=c(15,40)){

  seeds.id <- length(simpact.trans.net)

  # Add age at sampling
  new.transm.tab <- vector("list", length(seeds.id))

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

  ## Binding together all selected transmission transmission networks ##

  if(length(IDs.transm)>=1){

    for (q in 1:length(IDs.transm)){

      p <- IDs.transm[q]

      if(q==1){
        trans.sum <- new.transm.tab[[p]]
        rename.id <- paste0(p,".",trans.sum$id,".C")
        trans.sum$id <- rename.id
        trans.sum.rename.id <- trans.sum
      }
      else{

        read.trans.sum <- new.transm.tab[[p]]
        rename.id.read <- paste0(p,".",read.trans.sum$id,".C")
        read.trans.sum$id <- rename.id.read
        trans.sum.rename.id <- rbind(trans.sum.rename.id, read.trans.sum)
      }
    }


    trans.sum.men <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$GenderRec=="0")

    trans.sum.women <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$GenderRec=="1")

    trans.sum.age.limit.men <- dplyr::filter(trans.sum.men, trans.sum.men$age.i>=age.men[1], trans.sum.men$age.i<=age.men[2])

    trans.sum.age.limit.women <- dplyr::filter(trans.sum.women, trans.sum.women$age.i>=age.women[1], trans.sum.women$age.i<=age.women[2])

    perc.100 <- nrow(trans.sum.age.limit.men) + nrow(trans.sum.age.limit.women) # total number of individuals with age limit

    perc.seq.coverage <- round(perc.100*seq.cov/100) # total number of wanted individuals at seq.cov sequence coverage

    nrow.men <- nrow(trans.sum.age.limit.men) # number of available men

    nrow.women <- nrow(trans.sum.age.limit.women) # number of available women

    men.seq.coverage <- round((perc.seq.coverage*perc.men)/100) # number of wanted men

    women.seq.coverage <- perc.seq.coverage - men.seq.coverage # numbe rof wanted women

    diff.men <- nrow.men - men.seq.coverage # difference between number of available and wanted men

    diff.women <- nrow.women - women.seq.coverage # difference between number of available and wanted women

    ratio.emp <- perc.men/(100-perc.men)

    if(diff.men >0 & diff.women >0){ # perfect case

      samp.men <- sample(trans.sum.age.limit.men$id, men.seq.coverage)
      samp.women <- sample(trans.sum.age.limit.women$id, women.seq.coverage)
      ratio.seq <- men.seq.coverage/women.seq.coverage

    } else if(diff.men <0 & diff.women >0){ # we have less men than what we need

      samp.men <- sample(trans.sum.age.limit.men$id, nrow.men) # take all we have
      samp.women <- sample(trans.sum.age.limit.women$id, women.seq.coverage)
      ratio.seq <- nrow.men/women.seq.coverage

    } else if(diff.men >0 & diff.women <0){ # we have less women than what we need

      samp.men <- sample(trans.sum.age.limit.men$id, men.seq.coverage)
      samp.women <- sample(trans.sum.age.limit.women$id, nrow.women) # take all we have
      ratio.seq <- men.seq.coverage/nrow.women

    } else {   #if(diff.men <0 & diff.women <0){ # we have less men and women than what we need
      samp.men <- sample(trans.sum.age.limit.men$id, nrow.men) # take all we have
      samp.women <- sample(trans.sum.age.limit.women$id, nrow.women) # take all we have
      ratio.seq <- nrow.men/nrow.women

    }


    outputvector <- c(samp.men, samp.women)

    outputvec.stat <- list()

    outputvec.stat$outputvector <- outputvector
    outputvec.stat$men <- men.seq.coverage
    outputvec.stat$women <- women.seq.coverage
    outputvec.stat$ratio.seq <- ratio.seq
    outputvec.stat$ratio.emp <- ratio.emp

  }else{


    outputvec.stat <- list()

    outputvec.stat$outputvector <- outputvector <- NA
    outputvec.stat$men <- men.seq.coverage <- NA
    outputvec.stat$women <- women.seq.coverage <- NA
    outputvec.stat$ratio.seq <- ratio.seq <- NA
    outputvec.stat$ratio.emp <- ratio.emp <- NA

  }

  return(outputvec.stat)

}
