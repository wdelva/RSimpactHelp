load("testrun.RData")

library(data.table)
library(shape)
library(reshape2)
library(nlme)
library(phylobase)

#library(sna)
library(ggplot2)
library(stringr)
library(plyr)

degree.df.maker <- function(datalist = datalist, survey.time = 10, window.width = 1, only.new = TRUE){
  # For the people in the datalist,
  # make a dataframe that has the variables of the ptable,
  ptable.df = datalist$ptable
  rtable.df = datalist$rtable
  # but in addition, also the number of relationships that were ongoing (if only.new == FALSE)
  # or that were newly formed (if only.new == TRUE) at some time in the survey window.
  # (E.g. The one-year period from 9 to 10 years into the simulation).
  limit.val<-survey.time-window.width
  num.occurrency<-0
  {if (only.new){
    rtable.df$New[rtable.df$FormTime>limit.val & rtable.df$DisTime>survey.time]<-1
    new.list <- table(rtable.df$New)
    num.occurrency<-new.list[names(new.list)==1]
  }
  else{
    rtable.df$Ongoing[rtable.df$FormTime<=limit.val & rtable.df$DisTime>survey.time]<-1
    ongoing.list <- table(rtable.df$Ongoing)
    num.occurrency<-ongoing.list[names(ongoing.list)==1]
  }
  } 
  
  #create a new id for dataframe ptable, concatenation of ID and Gender
  ptable.df$ID.Gender <- paste(ptable.df $ID, ptable.df $Gender, sep=".")
  setkey(ptable.df , ID.Gender)
  #Get the degree of man from rtable with two variables: Var1 which is the factor from 1 to 8130 
  #since there are 8130 degrees, and Frequencies contain the values of degrees
  m.degree <- data.table(data.frame(table(rtable.df$IDm)))
  #Add a key to dataframe m.degree, to finally have 3 columns
  m.degree$ID.Gender <- paste(m.degree$Var1, 0, sep=".")
  #Change columns names
  setnames(m.degree, c("IDm", "degree", "ID.Gender"))
  setkey(m.degree, ID.Gender)
  
  #Get the degree of woman from rtable with two variables: Var1 which is the factor from 1 to 8130 
  #since there are 8130 degrees, and Freq contains the values of degrees
  w.degree <- data.table(data.frame(table(rtable.df$IDw)))
  w.degree$ID.Gender <- paste(w.degree$Var1, 1, sep=".")
  setnames(w.degree, c("IDw", "degree", "ID.Gender"))
  setkey(w.degree, ID.Gender)
  
  #Merge datasets
  degree.df <- merge(ptable.df, m.degree, all.x=TRUE)
  w.degree.df <- merge(ptable.df, w.degree, all.x=TRUE)
  #
  degree.df$degree[degree.df$Gender==1] <-  w.degree.df$degree[w.degree.df$Gender==1]
  #Change the values of degree from NA to 0
  degree.df$degree[is.na(degree.df$degree)] <- 0
  #degree.df$degree[degree.df$degree==0] <- 0.0001
  
  return(list(degree.df=degree.df,num.occurrency=num.occurrency))
}
degree.df <- degree.df.maker(datalist = datalist, survey.time = 10, window.width = 1, only.new = TRUE)