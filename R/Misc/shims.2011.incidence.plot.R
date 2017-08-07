library(tidyr, ggplot2)

#Swaziland SHIMS 2011 incidence by gender
xmen <- c(0.8,6.6,21.3,36.6,47.0,45.5,42.5)
xwom <- c(14.3, 31.5, 46.7, 53.8, 49.1, 39.7, 31.6)
time <- 1:7

data.plot <- data.frame(cbind(time, xwom,xmen))
names(data.plot) <- c("time", "Women", "Men")
## wide to long format
data.plot <- data.plot %>%
  tidyr::gather(Gender, "Incidence", 2:3)

ggplot(data.plot,
       aes(x=time, y=Incidence, group=Gender, label = Incidence)) +
  geom_point() + geom_line() + aes(colour = Gender)  +
  xlim("18-19","20-24","25-29","30-34","35-39","40-44","45-49") +
  ylim(0,60) +
  geom_text(aes(label=Incidence), size=4, vjust=-.8, show.legend = FALSE) +
  labs(title = "2011 HIV Prevalence by Age Group in Swaziland",
       x = "Age Group", y = " HIV Prevalence") +  theme_bw() +
  theme(axis.text.x  = element_text(vjust=0.5, size=14),
        axis.title.x = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=14),
        axis.title.y = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5, size=16))


