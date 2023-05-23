
##################################################################
## main run for one change point (testing)
##################################################################
rm(list=ls())
ptm = proc.time()
setwd("E:\\R_code\\Das_BentReg\\")

## for testing
kapa <- 1.01   ## parameter in sSIC
M <- 1000      ## the number of randem intervals
Kmax <- 20     ## the max number of change points

## for data
beta00 <- c(1, 1.5, -3, 2) 
t00 = 0
modtype <- 2   ##  modtype = 1, IID case,  modetype = 2, heteroscedasticity case
etype <- 1     ## error type: etype = 1 for N(0,1),  etype = 2 for t_3,   etype = 3 for 0.9N(0,1) + 0.1t_3, 
               ##  etype = 3 for Cauchy(0,1), etype = 5 for 0.9N(0,1) + 0.1Cauchy(0,1)
n = 200        ## sample size
NS <- 1000        ## simulation times

TAU <- 0.2 #seq(0.1, 0.2, by=0.1)

result <- matrix(NA, length(TAU), 3)
all.out <- matrix(NA, NS, length(TAU))
for(k in 1:length(TAU)){
  tau <- TAU[k]
  source("test_one_wbs.R")
  all.out[, k] <- est.number
  
  result[k,] <- c(mean(est.number, na.rm=TRUE),length(est.number[est.number==1])/length(subset(est.number,est.number>=0)),sd(est.number, na.rm=TRUE))
}
result
folder="E:\\paper_simulation\\Das_BentReg\\bent line quantile\\AllNewRsults\\"
filename = paste("test1","etype", etype, "modtype", modtype,"tau",TAU,"n", n, "NS",sep="-") 
write.csv(all.out, file=paste(folder, filename, ".csv", sep=''),  row.names=T)

colnames(result) <- c("mean","ratio","sd")
filename = paste("test1","etype", etype, "modtype", modtype, "tau",TAU,"n", n, sep="-") 
write.csv(result, file=paste(folder, filename, ".csv", sep=''),  row.names=T)
proc.time()-ptm




######################################################################
## main run for one change point (estimate)
######################################################################

rm(list=ls())
ptm = proc.time()
setwd("E:\\R_code\\Das_BentReg\\")

beta00 <- c(1, 1.5, -3, 2) 
t00 = 0
n = 200       # sample size
NS <- 1000
modtype <- 2  ##  modtype = 1, IID case,      modetype = 2, heteroscedasticity case
etype <- 3    ##  etype = 1 for N(0,1),       etype = 2 for t_3,        etype = 3 for 0.9N(0,1) + 0.1t_3, 
              ##  etype = 4 for Cauchy(0,1),  etype = 5 for 0.9N(0,1) + 0.1Cauchy(0,1)

TAU <- c(0.1,0.2,0.4,0.6,0.8,0.9) ##seq(0.6, 0.9, by=0.1)

for(k in 1:length(TAU)){
  tau <- TAU[k]
  source("one_changepoint.R")
}

proc.time()-ptm









