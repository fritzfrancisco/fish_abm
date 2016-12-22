#setwd("~/rolling_int_uptake_speed5")

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

r_data_all = read.csv(paste("rolling_int_uptake_speed5_",1,".csv",sep=""),header = FALSE)
colnames(r_data_all)=paste(1,colnames(r_data_all),sep="")

for(i in 2:15){
  r_data_i = read.csv(paste("rolling_int_uptake_speed5_",i,".csv",sep=""),header = FALSE)
  colnames(r_data_i)=paste(i,colnames(r_data_i),sep="")
  r_data_all = cbind(r_data_all,r_data_i)
}

nr_data_all = read.csv(paste("norolling_int_uptake_speed5_",1,".csv",sep=""),header = FALSE)
colnames(nr_data_all)=paste(1,colnames(nr_data_all),sep="")

for(i in 2:15){
  nr_data_i = read.csv(paste("norolling_int_uptake_speed5_",i,".csv",sep=""),header = FALSE)
  colnames(nr_data_i)=paste(i,colnames(nr_data_i),sep="")
  nr_data_all = cbind(nr_data_all,nr_data_i)
}

hist(conc,breaks = seq(0,80,by=4),xlim=c(0,70),main="Non-Rolling",xlab="cumulative intake")
hist(concI,breaks = seq(0,80,by=4),xlim=c(0,70),main="Rolling",xlab="cumulative intake")

list <- list(c(nr_data_all[500,]),c(r_data_all[500,]))
conc <- c(as.numeric(nr_data_all[500,]))
concI <- c(as.numeric(r_data_all[500,]))
var.test(conc,concI)

