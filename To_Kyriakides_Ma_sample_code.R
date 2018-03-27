#################################################
## 1) Data Manipulation: extracting NYC taxi data
## Description: The data includes variables:
##              drop-off time, pick-up time
##              longitude, latitude
#################################################
taxi <- read.table("data/NYC_yellow_2015_7_clean_small.csv", header = TRUE)
hist(taxi$pickup_latitude, breaks = 10)
plot(density(taxi$pickup_latitude))
hist(taxi$pickup_longitude)
plot(density(taxi$pickup_longitude))

taxi$tpep_pickup_datetime <- as.character(taxi$tpep_pickup_datetime)
taxi$tpep_dropoff_datetime <- as.character(taxi$tpep_dropoff_datetime)
dtparts <- t(as.data.frame(strsplit(taxi$tpep_pickup_datetime, ' ')))
row.names(dtparts) = NULL
thetimes = chron(dates=dtparts[,1],times=dtparts[,2], format=c('y-m-d','h:m:s'))
head(as.POSIXct(thetimes))
head(as.POSIXct(taxi$tpep_dropoff_datetime))

zones <-readOGR("/Users/Colin/Desktop/taxi_zones", layer="taxi_zones")
ggplot() +  geom_polygon(data=zones, aes(x=long, y=lat, group=group))



############################################
## 2) Simulation: multivariate Gaussian MCMC
############################################
MH <- function (Nsim, sigma, theta.init){
  theta <- theta.init
  samples <- c()
  ThetaNew <- function(theta, sigma){
    #cand: candidate point
    #ratio: acceptance ratio
    cand <- rnorm(1, mean=theta, sd= sigma)
    ratio <- (dnorm(cand, 2, 1)*dnorm(theta, mean=cand, sigma))/(dnorm(theta, 2, 1)*dnorm(cand, mean=theta, sigma))
    if (runif(1) <= ratio) cand
    else theta
  }
  for (i in 1:Nsim){
    theta <- ThetaNew(theta=theta, sigma=sigma)
    samples[i] <- theta
  }
  return(samples)
}
#Metropolis-Hastings J=N(theta', 1)
x <- MH(Nsim=5000, sigma=1, theta.init=2)
par(mfrow=c(3,1))
plot(x, type="l")
hist(x, main="histogram p(x)=N(2,1), J=N(theta', 1)")
plot(cumsum(x^3)/1:5000, type="l", xlab="Iterations", ylab=" ")
abline(a=14, b=0, col="red")
#integration estimate
mean(x^3)


######################################
## 3) LASSO: 10-fold cross validation
######################################
library(glmnet)
dat <- read.csv("nj.csv", header=TRUE)
delete <- c("serialno", "rectype", "STATE", "region", "division", "puma5", "PUMA1", "UNITTYPE", "VACSTAT", "TENURE", "BUSINES", "RENT", "PNUM", "RELATE", "GRADE", "OCCCEN5", "OCCSOC5")
dat <- dat[-match(delete, names(dat))] #remove columns
dat <- dat[dat$FINC>0,] # remove negative FINC
median(dat$FINC[dat$PERSONS==4])*0.8  #cut-off=60,000
inct = rep(60000, length(cut))
inct[cut<0] = (1+0.1*cut[cut<0])*60000
inct[cut>0] = (1+0.08*cut[cut>0])*60000
#fill in NA's with zero
dat$MORTG1[which(is.na(dat$MORTG1))] <-0
which(is.na(dat$MORTG1))
dat$MRT1AMT[which(is.na(dat$MRT1AMT))] <-0
dat$MORTG2[which(is.na(dat$MORTG2))] <-0

n <- nrow(dat)
fold <- sample(1:10, n, replace=TRUE)
xx <- as.matrix(dat)
y <- lowinc
yy = y*0
err <- c()
for(i in 1:10) {
  index <- which(fold!=i)
  hh=predict(glmnet(xx[-index,],alpha=0.5,factor(y)[-index],family="binomial"),newx=xx[,]) 
  uu=t(apply(hh[-index,],2,frate2,y[-index]))
  kk=which.min(uu[,2])
  yy[index]=hh[index,kk]>=uu[kk,1]
  err[i] <- table(y,yy)[1,2]+table(y,yy)[2,1]
}
table(y,yy)

cv.perf <- mean(err/108246)
