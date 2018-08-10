#' A function that returns age mixing patterns quantities in transmission networks
#' in scenarios where individuals are missing at completly at random
#' @param simpact.trans.net a list of transmission networks produced by \code{\link{transm.network.builder}}
#' @param limitTransmEvents Number of minimum transmission events to be considered in each transmission networks
#' @param timewindow Time interval
#' @param seq.cov Percentage of individulas considered for this transmission pattern scenario
#' @param age.group.15.25 age group between 15 and 25 years old
#' @param age.group.25.40 age group between 25 and 40 years old
#' @param age.group.40.50 age group between 40 and 50 years old
#' @return a vector of number of men and women in different age group, number of transmissions within all age groups, and mean and SD of age different between infectors and infectees
#' @examples
#' w <- CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#'                            limitTransmEvents = 7,
#'                            timewindow = c(30,40),
#'                            seq.cov = 70,
#'                            age.group.15.25 = c(15,25),
#'                            age.group.25.40 = c(25,40),
#'                            age.group.40.50 = c(40,50))

#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
#'


# true we record of infection - infection time

CAR.groups.fun.agemix <- function(simpact.trans.net = simpact.trans.net,
                                  limitTransmEvents = 7,
                                  timewindow = c(30,40),
                                  seq.cov = 70,
                                  age.group.15.25 = c(15,25),
                                  age.group.25.40 = c(25,40),
                                  age.group.40.50 = c(40,50)){


  seeds.id <- length(simpact.trans.net)

  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)

  for(i in 1:seeds.id){

    transm.age.df.ic <- as.data.frame(simpact.trans.net[[i]])

    age.samp.Rec <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBRec
    age.samp.Don <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBDon

    transm.age.i <- cbind(transm.age.df.ic, age.samp.Rec, age.samp.Don)

    new.transm.tab[[i]] <- transm.age.i

  }

  # id of people who got infection by seed event: seeds.id
  trans.network <- new.transm.tab

  seeds.id <- length(trans.network)


  ID.select <- vector() # ID of selected transmission network
  ID.select.count <- vector() # number of individuals in these networks

  for (i in 1: seeds.id) {


    trans.network.i <- as.data.frame(trans.network[[i]])

    if(nrow(trans.network.i)>=limitTransmEvents){


      ID.select <- c(ID.select, i)
      ID.select.count <- c(ID.select.count, nrow(trans.network.i))

    } # X if

  } # Y for


  infectionTable <- vector("list", length(ID.select))

  for(j in 1:length(ID.select)){

    p <- ID.select[j]

    trans.network.i <- as.data.frame(trans.network[[p]])

    trans.network.i <- trans.network.i[-1,]

    id.lab <- paste0(p,".",trans.network.i$id,".C")

    trans.network.i$id.lab <- id.lab

    infectionTable[[p]] <- trans.network.i
  }


  infecttable <- rbindlist(infectionTable)


  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow,
                             seq.cov = seq.cov,
                             age.limit = age.group.40.50[2])



  data.transm.agemix <- dplyr::filter(infecttable, infecttable$id.lab%in%mCAr.IDs)


  sort.partners.fun <- function(partner.table = partner.table){ # for receivers

    # age and gender structured receiver individuals

    num.15.25.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                              error=function(e) return(NULL))

    num.15.25.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                                error=function(e) return(NULL))


    num.25.40.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                              error=function(e) return(NULL))

    num.25.40.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                                error=function(e) return(NULL))


    num.40.50.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                              error=function(e) return(NULL))

    num.40.50.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                                error=function(e) return(NULL))

    # consider filter == men

    part.men.15.25.women.15.25 <- tryCatch(dplyr::filter(num.15.25.men, num.15.25.men$age.samp.Don >= age.group.15.25[1] & num.15.25.men$age.samp.Don < age.group.15.25[2]),
                                           error=function(e) return(NULL)) # table of women partners of men between 15 and 25 years old

    part.men.15.25.women.25.40 <- tryCatch(dplyr::filter(num.15.25.men, num.15.25.men$age.samp.Don >= age.group.25.40[1] & num.15.25.men$age.samp.Don < age.group.25.40[2]),
                                           error=function(e) return(NULL))

    part.men.15.25.women.40.50 <- tryCatch(dplyr::filter(num.15.25.men, num.15.25.men$age.samp.Don >= age.group.40.50[1] & num.15.25.men$age.samp.Don < age.group.40.50[2]),
                                           error=function(e) return(NULL))


    part.men.25.40.women.15.25 <- tryCatch(dplyr::filter(num.25.40.men, num.25.40.men$age.samp.Don >= age.group.15.25[1] & num.25.40.men$age.samp.Don < age.group.15.25[2]),
                                           error=function(e) return(NULL))

    part.men.25.40.women.25.40 <- tryCatch(dplyr::filter(num.25.40.men, num.25.40.men$age.samp.Don >= age.group.25.40[1] & num.25.40.men$age.samp.Don < age.group.25.40[2]),
                                           error=function(e) return(NULL))

    part.men.25.40.women.40.50 <- tryCatch(dplyr::filter(num.25.40.men, num.25.40.men$age.samp.Don >= age.group.40.50[1] & num.25.40.men$age.samp.Don < age.group.40.50[2]),
                                           error=function(e) return(NULL))


    part.men.40.50.women.15.25 <- tryCatch(dplyr::filter(num.40.50.men, num.40.50.men$age.samp.Don >= age.group.15.25[1] & num.40.50.men$age.samp.Don < age.group.15.25[2]),
                                           error=function(e) return(NULL))

    part.men.40.50.women.25.40 <- tryCatch(dplyr::filter(num.40.50.men, num.40.50.men$age.samp.Don >= age.group.25.40[1] & num.40.50.men$age.samp.Don < age.group.25.40[2]),
                                           error=function(e) return(NULL))

    part.men.40.50.women.40.50 <- tryCatch(dplyr::filter(num.40.50.men, num.40.50.men$age.samp.Don >= age.group.40.50[1] & num.40.50.men$age.samp.Don < age.group.40.50[2]),
                                           error=function(e) return(NULL))


    N.partners <- c(nrow(part.men.15.25.women.15.25), nrow(part.men.15.25.women.25.40), nrow(part.men.15.25.women.40.50),
                    nrow(part.men.25.40.women.15.25), nrow(part.men.25.40.women.25.40), nrow(part.men.25.40.women.40.50),
                    nrow(part.men.40.50.women.15.25), nrow(part.men.40.50.women.25.40), nrow(part.men.40.50.women.40.50))

    return(N.partners)

  }


  trans.sum.age.limit <- data.transm.agemix

  # Group 15 - 25
  ###############

  trans.sum.men.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.samp.Rec >= age.group.15.25[1] & trans.sum.age.limit$age.samp.Rec < age.group.15.25[2])

  trans.sum.women.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.samp.Rec >= age.group.15.25[1] & trans.sum.age.limit$age.samp.Rec < age.group.15.25[2])


  # Group 25 - 40
  ###############

  trans.sum.men.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.samp.Rec >= age.group.25.40[1] & trans.sum.age.limit$age.samp.Rec < age.group.25.40[2])

  trans.sum.women.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.samp.Rec >= age.group.25.40[1] & trans.sum.age.limit$age.samp.Rec < age.group.25.40[2])



  # Group 40 - 50
  ###############

  trans.sum.men.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.samp.Rec >= age.group.40.50[1] & trans.sum.age.limit$age.samp.Rec < age.group.40.50[2])

  trans.sum.women.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.samp.Rec >= age.group.40.50[1] & trans.sum.age.limit$age.samp.Rec < age.group.40.50[2])



  partners.age.str <- sort.partners.fun(data.transm.agemix) # 154

  ouput.transm.dat <- c(nrow(trans.sum.men.15.25), nrow(trans.sum.women.15.25),
                        nrow(trans.sum.men.25.40), nrow(trans.sum.women.25.40),
                        nrow(trans.sum.men.40.50), nrow(trans.sum.women.40.50),

                        partners.age.str)


  # Age difference statistics #
  #############################
  AD <- abs(abs(data.transm.agemix$TOBDon) - abs(data.transm.agemix$TOBRec))
  mean.AD <- mean(AD)
  med.AD <- median(AD)
  sd.AD <- sd(AD)

  # Mixed effect models #
  #######################
  # fit.agemix.trans.women <- fit.agemix.trans.women(datatable = data.transm.agemix)
  # fit.agemix.trans.men <- fit.agemix.trans.men(datatable = data.transm.agemix)


  ouput.transm.dat.AD <- c(ouput.transm.dat, mean.AD, med.AD, sd.AD)


  val.names <- c("num.men.15.25", "num.women.15.25",
                 "num.men.25.40", "num.women.25.40",
                 "num.men.40.50", "num.women.40.50",

                 "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                 "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                 "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",

                 "mean.AD", "median.AD", "sd.AD")


  names(ouput.transm.dat.AD) <- val.names


  return(ouput.transm.dat.AD)
}

### ----------------------------------------

#
# w <- CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#                            limitTransmEvents = 7,
#                            timewindow = c(30,40),
#                            seq.cov = 70,
#                            #    seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
#                            age.group.15.25 = c(15,25),
#                            age.group.25.40 = c(25,40),
#                            age.group.40.50 = c(40,50))




