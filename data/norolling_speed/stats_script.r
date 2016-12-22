#setwd("~/norolling_speed")

get_rates <- function(data_set){
  rates <- NULL
  start <- NULL
  for(i in 1:ncol(data_set)){
    x = seq(0, nrow(data_set) - 1, by = 1)
    y = data_set[,i]
    mod = lm(y ~ x)
    rates[i] <- mod[[1]][2]
    start[i] <- mod[[1]][1]
  }
  return(data.frame(rates,start))
}

nrspeed0 <- data.frame(t(read.csv("norolling_speed2.5.csv", header = FALSE)))
row.names(nrspeed0) <- NULL

nrspeed1 <- data.frame(t(read.csv("norolling_speed5.csv", header = FALSE)))
row.names(nrspeed1) <- NULL

nrspeed2 <- data.frame(t(read.csv("norolling_speed7.5.csv", header = FALSE)))
row.names(nrspeed2) <- NULL

nrspeed3 <- data.frame(t(read.csv("norolling_speed10.csv", header = FALSE)))
row.names(nrspeed3) <- NULL

nrspeed4 <- data.frame(t(read.csv("norolling_speed12.5.csv", header = FALSE)))
row.names(nrspeed4) <- NULL

nrspeed5 <- data.frame(t(read.csv("norolling_speed15.csv", header = FALSE)))
row.names(nrspeed5) <- NULL

nrspeed6 <- data.frame(t(read.csv("norolling_speed17.5.csv", header = FALSE)))
row.names(nrspeed6) <- NULL

nrspeed7 <- data.frame(t(read.csv("norolling_speed20.csv", header = FALSE)))
row.names(nrspeed7) <- NULL

matplot(nrspeed1,type="l",xlab = "iterations",ylab = "environmental quality",main="Depletion: Non-Rolling",ylim=c(1500,3700))

rates_025<- get_rates(nrspeed0)
rates_050<- get_rates(nrspeed1)
rates_075<- get_rates(nrspeed2)
rates_100<- get_rates(nrspeed3)
rates_125<- get_rates(nrspeed4)
rates_150<- get_rates(nrspeed5)
rates_175<- get_rates(nrspeed6)
rates_200<- get_rates(nrspeed7)

boxplot(rates_025$rates,rates_050$rates,rates_075$rates,rates_100$rates,rates_125$rates,rates_150$rates,rates_175$rates,rates_200$rates,xlab="speed",ylab="rate of depletion")

shapiro.test(rates_025$rates)
shapiro.test(rates_050$rates)
shapiro.test(rates_075$rates)
shapiro.test(rates_100$rates)
shapiro.test(rates_125$rates)
shapiro.test(rates_150$rates)
shapiro.test(rates_175$rates)
shapiro.test(rates_200$rates)

nr_rate_list <- c(rates_025$rates,rates_050$rates,rates_075$rates,rates_100$rates,rates_125$rates,rates_150$rates,rates_175$rates,rates_200$rates)
nr_speed_list <- rep(seq(2.5,20,by=2.5),each=30)

boxplot(nr_rate_list ~ nr_speed_list,xlab="speed",ylab="depletion rate",ylim=c(-1,-0),main="Without Rolling Behaviour")
bartlett.test(rate_list~speed_list)
kruskal.test(rate_list~speed_list)
