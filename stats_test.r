library("ggplot2")

setwd("/home/fritz/Documents/Projects/Couzin Lab Projects/stats")


qualities = {}
qualities1 = {}
qualities2 = {}

times ={}
for(i in 1:10000){
  u = i %% 10
  if(u == 0){
    u = 10
  }
    t = trunc((i-1)/10) + 1
    
    qualities[i] = sq[t,u]
    qualities1[i] = sq2[t,u]
    qualities2[i] = sq3[t,u]
    
    times[i] = trunc(i/10)+1
}

trial = data.frame(times,qualities,qualities1,qualities2)

ggplot(data = trial) + 
  geom_smooth(mapping = aes(x = trial$times,y = trial$qualities),method = "lm") +
  geom_smooth(mapping = aes(x = trial$times,y = trial$qualities1),method = "lm") +
  geom_smooth(mapping = aes(x = trial$times,y = trial$qualities2),method = "lm")
 # geom_point(size = 0.05)

model <- lm(qualities ~ times)
summary(model)

# feeding rate of single fish per minute
model[[1]][2] * 2000/45
