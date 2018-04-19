
#remove columns
#Using MPI on the cluster
#mpirun -np 4 Rscript R/Misc/Pre.hhohhoMaxARTFinal.Sim/T.pre.hho.simpact.wrapper.easyABC.run.R


#df <- subset(df, select = -c(x,y))
#drop = c("x","y")
#df <- df[,!(names(df) %in% drop)]

#new.xdesign.ssd.pc.c <- new.xdesign.ssd.pc.c[order(new.xdesign.ssd.pc.c$ssd.c),]
#inputANDoutput.select <- inputANDoutput.select %>% dplyr::group_by(sim.id) %>% dplyr::summarise_each(funs(mean))

#pairs(inputANDoutput.select[, x.variables[10:14]], col = 1+inputANDoutput.select$is.complete, pch = 16, cex = 2)#,

#rename col
#names(df)[names(df) == 'old.var.name'] <- 'new.var.name'

#addline_format <- function(x,...){
#  gsub('\\s','\n',x)
#}

#qplot(ICD1_TB_General, data = plotme_now, fill = hiv_status, ylab = "Frequency", xlab = "Top 10 causes of death with HIV status") +
#  scale_x_discrete(labels=addline_format(c("A_2" = "Tuberculosis","J_5" = "Pneumonia","I_25" = "Stroke","N_6" = "Kidney Failure","A_8" = "Other Sepsis","G_5" = "Meningitis", "B_4" = "HIV Infection","A_1" = "Intestinal Infectious" ,"I_1" = "Hypertension" ,"B_9" = "Crypto coccosis"))) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#visualise the correlation matrix
#The asterisks indicate the significance levels of the correlations.
#Each significance level is associated to a symbol :
#p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
#library(PerformanceAnalytics)
#chart.Correlation(mydata, histogram=TRUE, pch=19)

#rle produces two vectors (lengths and values ). The length of values vector gives you the number of unique values.
#x<-c(1,2,3,1,2,3,4,6)
#length(rle(sort(x))$values)

#values in a list
#vc <- c(1,2,3,10)
#dt[dt$id %in% vc,]
#dplyr::filter(df, fct %in% vc)


#Sort data.frame
#dd <- dd[order(dd$b, decreasing = FALSE),] one variable
#library(dplyr)
#arrange(mtcars, mpg, cyl, wt) by more than one column

#Get the first and last row from grouped data.
#df %>%
#  group_by(id) %>%
#  arrange(stopSequence) %>% #you can arrange by other columns
#  filter(row_number()==1 | row_number()==n())

#Create new index using the grouped data
#df %>%
#  group_by(id) %>%
#  mutate(range = 1:n())

#get summary on grouped data
#df %>% group_by(group) %>% summarize(mean=mean(dt), sum=sum(dt))

#from psych
#describeBy(df$dt, df$group, mat = TRUE)

#On the terminal
#history | grep YOUR_STRING

#Check if column contains something
#df <- iris[grep("what", df$colmn), ]


#convert the df in to a vector with no df names
#unname(unlist(df[1,]))

#get a vector
#unique(unlist(x, use.names = FALSE))

#round off
#sprintf("%02d",file.read)



# inc.plot.t.vm <- subset(inci.select.sum.plot, Gender=="Men")
# names(inc.plot.t.vm)[names(inc.plot.t.vm) == 'Value'] <- 'inc.m'
#
# inc.plot.t.vw <- subset(inci.select.sum.plot, Gender=="Women")
# names(inc.plot.t.vw)[names(inc.plot.t.vw) == 'Value'] <- 'inc.w'
#
# prev.plot.t.vm <- subset(prev.select.sum.plot, Gender=="Men")
# names(prev.plot.t.vm)[names(prev.plot.t.vm) == 'Value'] <- 'prev.m'
#
# prev.plot.t.vw <- subset(prev.select.sum.plot, Gender=="Women")
# names(prev.plot.t.vw)[names(prev.plot.t.vw) == 'Value'] <- 'prev.w'
#
# data <- as.data.frame(cbind(prev.plot.t.vm$x.lab, inc.plot.t.vm$inc.m, inc.plot.t.vw$inc.w,
#                             prev.plot.t.vm$prev.m, prev.plot.t.vw$prev.w))
#
# names(data) <- c("age.group.point", "inc.m", "inc.w", "prev.m", "prev.w")
#
# data$inc.m[13:14] <- c(NA,NA)
#
# obj1 <- xyplot(inc.m + inc.w  ~ age.group.point, data,  type = c("l","p"), pch=19, lwd=2,
#                col =c("red","blue"), ylab = "HIV Incidence", panel=function(x, y,...) {
#                  panel.xyplot(x, y,...)
#                  panel.text(data$age.group.point, data$inc.w,
#                             labels=round(data$inc.w, digits=2), pos=2, offset=.8)
#                  panel.text(data$age.group.point, data$inc.m,
#                             labels=round(data$inc.m, digits=2), pos=4, offset=.8)
#                })
#
# obj2 <- xyplot(prev.m + prev.w ~ age.group.point, data, type = c("l","p"), pch=19, lwd=2,
#                col =c("green","purple"), ylab = "HIV Prevalence", panel=function(x, y,...) {
#                  panel.xyplot(x, y,...)
#                  panel.text(data$age.group.point, data$prev.w,
#                             labels=round(data$prev.w, digits=2), pos=4, offset=.8)
#                  panel.text(data$age.group.point, data$prev.m,
#                             labels=round(data$prev.m, digits=2), pos=4, offset=.8)
#                } )
#
# update(doubleYScale(obj1, obj2, style1 = 0, style2 = 1, add.ylab2 = TRUE,
#                     auto.key=list(title="", text = c("Incidence Women", "Incidence Men", "Prevalance Women", "Prevalence Man"),
#                                   columns = 4)),
#        par.settings = simpleTheme(col = c("blue", "red", "green","purple")),
#        scales = list(x=list(at=c(1,2, 3, 4, 5, 6, 7), labels = c("18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"))),
#        xlab= "Age Group", lty = 1:4)

#
# get a colored graph
# ## Add a new column that will decide this
# dfm$ismatch<-ifelse(
#   with(dfm,interaction(variable, value)) %in%
#     with(data, interaction(position,type)),
#   "match","nomatch")
# ### draw thsi graph with colors
# ggplot(dfm, aes(x=variable, y=gene, label=value, fill=ismatch)) +
#   geom_text(colour="black") +
#   geom_tile(alpha=0.5)

#zoomin function from plotrix
#zoomInPlot(rnorm(100),rnorm(100),rxlim=c(-0.5,0.5),rylim=c(-1,1),
           #            zoomtitle="Zoom In Plot",titlepos=-1.5)

#set date diff
#as.Date("2015-12-01")
#difftime(EndDate, BeginDate,  units = "days") #Weeks, Months


# Compare functions
# library(dplyr)
# library(microbenchmark)
#
# # Original example
# microbenchmark(
#   df1<-subset(airquality, Temp>80 & Month > 5),
#   df2<-filter(airquality, Temp>80 & Month > 5)
# )



###################Viral Load Plots ################### IGNORE
# yaxis <- vl.cutoff + 0.2 # add the axis to make visual
#
# #Visualise the VL points
# q <- ggplot()
# q <- q + geom_point(data=VLevent.df, aes(x=Time, y=Log10VL, group=ID, colour = Desc))
# #q <- q + geom_line(data=VLevent.df, aes(x=Time, y=Log10VL, group=ID, colour = Desc))
# q <- q + geom_hline(yintercept=vl.cutoff) +
#   annotate("text", datalist$itable$hivseed.time[1], vl.cutoff + 0.2, label = "  VLCutoff")
# q <- q + geom_vline(xintercept=timepoint-0.5, colour = "blue") +
#   annotate("text", timepoint, max(VLevent.df$Log10VL), label = "  TimePoint")
# q <- q +theme_bw() + theme(legend.position = "right") +
#   theme(panel.background = element_blank(),panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.y = element_blank(), panel.ontop = TRUE) +
#   theme(legend.background = element_rect(colour = "black"))

#Get the last recorded VL and the desc
#if StartedART and below vl.cutoff then suppressed otherwise NOT suppressed

##VLevent.df <- dplyr(VLevent.df,.(ID), tail,1)
##VLevent.df <- VLevent.df[, .SD[c(.N)], by=ID]
#####################################





