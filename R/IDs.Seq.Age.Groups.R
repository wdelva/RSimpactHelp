
#' Subset sequences based on gender ratio and age of individuals
#'
#' @param simpact.trans.net Transmission networks computed by \code{\link{transmission.network.builder()}}
#' @param limitTransmEvents Choose transmission network with at least this  amount of individuals
#' @param timewindow Time interval in which sampling are considered
#' @param seq.cov Sequece coverage
#' @param seq.gender.ratio Gender imbalance ratio (women/(women+men))
#' @param age.group.15.25 Age group
#' @param age.group.25.40 Age group
#' @param age.group.40.50 Age group
#' @return A vector of selected IDs
#' @importFrom dplyr filter
#' @export
#'


IDs.Seq.Age.Groups <- function(simpact.trans.net = simpact.trans.net,
                               limitTransmEvents = 7,
                               timewindow = c(10,40),
                               seq.cov = 70,
                               seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                               age.group.15.25 = c(15,25),
                               age.group.25.40 = c(25,40),
                               age.group.40.50 = c(40,50)){

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

    trans.sum.age.limit <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$age.i<=age.group.40.50[2])

    trans.sum.age.limit <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$SampTime >= timewindow[1] & trans.sum.age.limit$SampTime <= timewindow[2])

    # Group 15 - 25
    ###############

    trans.sum.men.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.i >= age.group.15.25[1] & trans.sum.age.limit$age.i < age.group.15.25[2])

    trans.sum.women.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.i >= age.group.15.25[1] & trans.sum.age.limit$age.i < age.group.15.25[2])

    perc.100.15.25 <- nrow(trans.sum.men.15.25) + nrow(trans.sum.women.15.25) # total number of individuals with age limit

    trans.sum.men.women.15.25 <- rbind(trans.sum.men.15.25, trans.sum.women.15.25)


    perc.women.15.25 <- round(perc.100.15.25 * seq.cov * seq.gender.ratio/100)

    perc.men.15.25 <- round(perc.100.15.25 * seq.cov * (1-seq.gender.ratio)/100)

    if(perc.men.15.25 <= length(trans.sum.men.15.25$id)){

      x.id.15.25 <- sample(trans.sum.men.15.25$id, perc.men.15.25)

      }else{

        x.id.15.25 <- sample(trans.sum.men.15.25$id, length(trans.sum.men.15.25$id))

      }

    if(perc.women.15.25 <= length(trans.sum.women.15.25$id)){

      y.id.15.25 <- sample(trans.sum.women.15.25$id, perc.women.15.25)

    }else{

      y.id.15.25 <- sample(trans.sum.women.15.25$id, length(trans.sum.women.15.25$id))

    }


    samp.all.15.25 <- c(x.id.15.25, y.id.15.25)




    # Group 25 - 40
    ###############

    trans.sum.men.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.i >= age.group.25.40[1] & trans.sum.age.limit$age.i < age.group.25.40[2])

    trans.sum.women.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.i >= age.group.25.40[1] & trans.sum.age.limit$age.i < age.group.25.40[2])

    perc.100.25.40 <- nrow(trans.sum.men.25.40) + nrow(trans.sum.women.25.40) # total number of individuals with age limit

    trans.sum.men.women.25.40 <- rbind(trans.sum.men.25.40, trans.sum.women.25.40)

    perc.women.25.40 <- round(perc.100.25.40 * seq.cov * seq.gender.ratio/100)

    perc.men.25.40 <- round(perc.100.25.40 * seq.cov * (1-seq.gender.ratio)/100)


    if(perc.men.25.40 <= length(trans.sum.men.25.40$id)){

      x.id.25.40 <- sample(trans.sum.men.25.40$id, perc.men.25.40)

    }else{

      x.id.25.40 <- sample(trans.sum.men.25.40$id, length(trans.sum.men.25.40$id))

    }

    if(perc.women.25.40 <= length(trans.sum.women.25.40$id)){

      y.id.25.40 <- sample(trans.sum.women.25.40$id, perc.women.25.40)

    }else{

      y.id.25.40 <- sample(trans.sum.women.25.40$id, length(trans.sum.women.25.40$id))

    }


    samp.all.25.40 <- c(x.id.25.40, y.id.25.40)


    # Group 40 - 50
    ###############

    trans.sum.men.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.i >= age.group.40.50[1] & trans.sum.age.limit$age.i < age.group.40.50[2])

    trans.sum.women.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.i >= age.group.40.50[1] & trans.sum.age.limit$age.i < age.group.40.50[2])

    perc.100.40.50 <- nrow(trans.sum.men.40.50) + nrow(trans.sum.women.40.50) # total number of individuals with age limit

    trans.sum.men.women.40.50 <- rbind(trans.sum.men.40.50, trans.sum.women.40.50)


    perc.women.40.50 <- round(perc.100.40.50 * seq.cov * seq.gender.ratio/100)

    perc.men.40.50 <- round(perc.100.40.50 * seq.cov * (1-seq.gender.ratio)/100)



    if(perc.men.40.50 <= length(trans.sum.men.40.50$id)){

      x.id.40.50 <- sample(trans.sum.men.40.50$id, perc.men.40.50)

    }else{

      x.id.40.50 <- sample(trans.sum.men.40.50$id, length(trans.sum.men.40.50$id))

    }

    if(perc.women.40.50 <= length(trans.sum.women.40.50$id)){

      y.id.40.50 <- sample(trans.sum.women.40.50$id, perc.women.40.50)

    }else{

      y.id.40.50 <- sample(trans.sum.women.40.50$id, length(trans.sum.women.40.50$id))

    }


    samp.all.40.50 <- c(x.id.40.50, y.id.40.50)


    samp.all <- c(samp.all.15.25, samp.all.25.40, samp.all.40.50)

  }else{
    samp.all <- NA
  }

  return(samp.all)
}
