#setwd("~/rolling_speed")

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

rspeed0 <- data.frame(t(read.csv("rolling_speed2.5.csv", header = FALSE)))
row.names(rspeed0) <- NULL

rspeed1 <- data.frame(t(read.csv("rolling_speed5.csv", header = FALSE)))
row.names(rspeed1) <- NULL

rspeed2 <- data.frame(t(read.csv("rolling_speed7.5.csv", header = FALSE)))
row.names(rspeed2) <- NULL

rspeed3 <- data.frame(t(read.csv("rolling_speed10.csv", header = FALSE)))
row.names(rspeed3) <- NULL

rspeed4 <- data.frame(t(read.csv("rolling_speed12.5.csv", header = FALSE)))
row.names(rspeed4) <- NULL

rspeed5 <- data.frame(t(read.csv("rolling_speed15.csv", header = FALSE)))
row.names(rspeed5) <- NULL

rspeed6 <- data.frame(t(read.csv("rolling_speed17.5.csv", header = FALSE)))
row.names(rspeed6) <- NULL

rspeed7 <- data.frame(t(read.csv("rolling_speed20.csv", header = FALSE)))
row.names(rspeed7) <- NULL

matplot(rspeed1,type="l",xlab = "iterations",ylab = "environmental quality",main="Depletion: Rolling",ylim=c(1500,3700))

rates_025<- get_rates(rspeed0)
rates_050<- get_rates(rspeed1)
rates_075<- get_rates(rspeed2)
rates_100<- get_rates(rspeed3)
rates_125<- get_rates(rspeed4)
rates_150<- get_rates(rspeed5)
rates_175<- get_rates(rspeed6)
rates_200<- get_rates(rspeed7)

boxplot(rates_025$rates,rates_050$rates,rates_075$rates,rates_100$rates,rates_125$rates,rates_150$rates,rates_175$rates,rates_200$rates,xlab="speed",ylab="rate of depletion")

shapiro.test(rates_025$rates)
shapiro.test(rates_050$rates)
shapiro.test(rates_075$rates)
shapiro.test(rates_100$rates)
shapiro.test(rates_125$rates)
shapiro.test(rates_150$rates)
shapiro.test(rates_175$rates)
shapiro.test(rates_200$rates)

r_rate_list <- c(rates_025$rates,rates_050$rates,rates_075$rates,rates_100$rates,rates_125$rates,rates_150$rates,rates_175$rates,rates_200$rates)
r_speed_list <- rep(seq(2.5,20,by=2.5),each=30)

boxplot(r_rate_list ~ r_speed_list,xlab="speed",ylab="depletion rate",ylim=c(-1,-0),main="Rolling behaviour")
bartlett.test(r_rate_list~r_speed_list)
kruskal.test(r_rate_list~r_speed_list)
