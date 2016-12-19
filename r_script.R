setwd("Documents/Studium/Master/VTKCouzinLab")

rolling <- data.frame(t(read.csv("rolling.csv", header = FALSE)))
row.names(rolling) <- NULL

no_rolling <- data.frame(t(read.csv("norolling.csv", header = FALSE)))
row.names(no_rolling) <- NULL

matplot(rolling, type = "l")
matplot(no_rolling,type = "l")

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

rolling_rates <- get_rates(rolling)
no_rolling_rates <- get_rates(no_rolling)

t.test(rolling_rates$rates, no_rolling_rates$rates)
boxplot(rolling_rates$rates, no_rolling_rates$rates)
