#' A function that recruits the MaxART clients
#'
#'
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata}}
#' @return the ptable datalist with clients tagged as study Yes/No
#'
#' @examples
#' data(datalist)
#' datalist.testme <- recruit.study.clients(datalist = datalist.test)
#' datalist
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

recruit.study.clients <- function(datalist = datalist){

  DT.recruit <- datalist$ptable
  recruit.time <- datalist$itable$maxart.starttime[1]


  DT.recruit.site <- DT.recruit %>%
    mutate(Age = recruit.time - TOB) %>%  #Assign age
    mutate(study.selected = "No") %>%     #Assign NOT selected to all clients
    filter(Age >= 18 &                    #greater than
          InfectTime !=Inf &              #HIV positive
          pfacility != "Not Hhohho" &     #study region
          !(ID %in% datalist$ttable$ID)   #treatment naive
          ) %>%
    as.data.frame

  #table(DT.recruit.site$pfacility)


 facility.name.prop <- list(c("Bulandzeni Clinic", 5), c("Herefords Clinic", 8),
                            c("Horo Clinic", 11), c("Maguga Clinic", 7),
                            c("Malandzela Nazarene Clinic", 3),
                            c("Mshingishingini Nazarene Clinic", 6),
                            c("Ndvwabangeni Nazarene Clinic", 8),
                            c("Pigg's Peak Government Hospital", 23),
                            c("Sigangeni Clinic", 6), c("Siphocosini Clinic", 7),
                            c("Hhukwini Clinic", 4), c("Ndzingeni Nazarene Clinic", 3),
                            c("Ntfonjeni Clinic", 7), c("Peak Nazarene Clinic", 2))



  #break the whole recruitment datalist for sampling by facility
  recruit.sample <- split(DT.recruit.site, DT.recruit.site$pfacility)

  for(i in 1:length(recruit.sample)){
    prop.name <- facility.name.prop[[i]][1]
    prop.prop <- as.numeric(facility.name.prop[[i]][2])/100

    choice.n <- round(prop.prop * length(recruit.sample[[prop.name]]$ID) ,0)

    recruit.sample[[prop.name]]$study.selected[which(recruit.sample[[prop.name]]$ID %in%
                                       sample(recruit.sample[[prop.name]]$ID, choice.n))] <- "Yes"


  }


  datalist.clients <- data.frame(do.call(rbind, recruit.sample))

  return(datalist.clients)
}

#datalist.testme <- recruit.study.clients(datalist.test)
#table(datalist.testme$study.selected)


#DT.recruit$study.selected[which(DT.recruit$ID %in% datalist.clients$ID)] <-
#  datalist.clients$study.selected[which(DT.recruit$ID %in% datalist.clients$ID)]



### Will need to make these target statistics


#Age = (female, total)
# recruit.age.bin <- list("A18.19" = c(91/100, 2/100), "A20.24" = c(84/100, 13/100),
#                         "A25.29" = c(65/100, 22/100), "A30.34" = c(57/100, 22/100),
#                         "A35.39" = c(53/100, 14/100), "A40.44" = c(59/100, 10/100),
#                         "A45.49" = c(68/100, 7/100), "A50.54" = c(65/100, 4/100),
#                         "A55.150" = c(53/100, 6/100))
#
#
# facility.age.aver.prop <- list("A18.19" = 2/100, "A20.24" = 13/100, "A25.29" = 22/100, "A30.34" = 20/100,
#                                "A35.39" = 15/100, "A40.44" = 10/100, "A45.49" = 7/100, "A50.54" = 5/100,
#                                "A55.150" = 6/100)

