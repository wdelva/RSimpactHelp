
#' Subset sequences based on gender ratio and age of individuals
#'
#' @param simpact.trans.net Transmission networks computed by \code{\link{transmission.network.builder()}}
#' @param limitTransmEvents Choose transmission network with at least this  amount of individuals
#' @param perc.men Precentage of men to consider when selecting the sequences
#' @param age.limit Age limit for all individuals
#' @return A list with a vector of ID's of selected sequences, number of men, and women for the selcted sequences and the ratio between men and women among the selcted sequences because in some settings we may have less or more female or male individuals
#' @importFrom dplyr filter
#' @export


select.indiv.seq.fun <- function(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 3,
                                 perc.men = 50, 
                                 age.limit=65){
  
  seeds.id <- length(simpact.trans.net)
  
  new.transm.tab <- vector("list", length(seeds.id))
  
  for(i in 1:seeds.id){
    
    transm.age.i <- as.data.frame(simpact.trans.net[[i]])
    
    age.i <- transm.age.i$SampTime - transm.age.i$TOBRec
    
    transm.age.i <- cbind(transm.age.i, age.i)
    
    new.transm.tab[[i]] <- transm.age.i
    
  }
  
  # Selceted networks
  
  IDs.transm <- vector()
  
  TransmEventsCountVector <- vector()
  
  for(k in 1:seeds.id){
    trans.net.i.check <- as.data.frame(new.transm.tab[[k]])
    
    if(nrow(trans.net.i.check)>=limitTransmEvents){
      
      TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))
      
      IDs.transm <- c(IDs.transm, k)
    }
  }
  
  ## Binding transmission transmission networks ##
  
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
  
  trans.sum.age.limit <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$age.i<=age.limit)
  
  trans.sum.men <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0")
  
  trans.sum.women <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1")  
  
  perc.100 <- nrow(trans.sum.men) + nrow(trans.sum.women) # 20 & 50
  
  nrow.men <- nrow(trans.sum.men) # number of available men
  
  nrow.women <- nrow(trans.sum.women) # number of available women
  
  men <- round((perc.100*per.men)/100) # number of wanted men
  
  women <- perc.100 - men # numbe rof wanted women
  
  diff.men <- nrow.men - men # difference between number of available and wanted men
  
  diff.women <- nrow.women - women # difference between number of available and wanted women
  
  ratio <- men/women
  
  if(diff.men >0 & diff.women >0){
    
    samp.men <- sample(trans.sum.men$id, men)
    samp.women <- sample(trans.sum.women$id, women)
    
  } else if(diff.men <0 & diff.women >0){
    
    samp.men <- sample(trans.sum.men$id, nrow.men)
    samp.women <- sample(trans.sum.women$id, women)
    
  } else if(diff.men >0 & diff.women <0){
    
    samp.men <- sample(trans.sum.men$id, men)
    samp.women <- sample(trans.sum.women$id, nrow.women) 
    
  } else {   #if(diff.men <0 & diff.women <0){
    samp.men <- sample(trans.sum.men$id, nrow.men)
    samp.women <- sample(trans.sum.women$id, nrow.women) 
    
  }
  
  outputvector <- c(samp.men, samp.women)
  
  outputvec.stat <- list()
  
  outputvec.stat$outputvector <- outputvector
  outputvec.stat$men <- men
  outputvec.stat$women <- women
  outputvec.stat$ratio <- ratio
  
  return(outputvec.stat)
  
}

