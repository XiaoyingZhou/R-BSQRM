##############################################################################
### simulation test code forthree thresholds linear quantile regression 
### test by WBS method,   
### used for paper simulation
### Author: Xiaoying Zhou  
### Date: On August 24, 2017
### 
###############################################################################
#rm(list=ls())
ptm <- proc.time()
#set.seed(2017)
library(quantreg)

## one change point------simulation data
bentdat <- function(n, beta00, t00, tau, modtype, etype){
  x <- runif(n, -4, 4)
  z <- rbinom(n, 1, 0.5)
  num.break <- length(t00)
  xx <- matrix(NA, n, num.break)
  for(i in 1:num.break){
    xx[,i] <- pmax(x-t00[i],0)
  }
  xz <- cbind(1, x, xx, z)
  if(etype ==1){
    e <- rnorm(n, 0, 0.5) - qnorm(tau, 0, 0.5)    # noraml distribution
  }else if(etype ==2){
    e <- rt(n, 4) - qt(tau, 4)            # t_3 distribution
  }else if(etype ==3){
    prop <- 0.1
    B <- rbinom(n, 1, prop)
    e <- (B == 0)*rnorm(n, 0, 0.5) + (B == 1)*rt(n, 4) -
      ((B == 0)*qnorm(tau, 0, 0.5) + (B == 1)*qt(tau, 4))    ## mixture 0.9*N(0,1) + 0.1 t_3
  }else if(etype ==4){
    e <- rcauchy(n, 0, 0.5) - qcauchy(tau, 0, 0.5)   ## Cauchy distribution
  } else if(etype ==5){
    prop <- 0.1
    B <- rbinom(n, 1, prop)
    e <-(B == 0)*rnorm(n, 0, 1) + (B == 1)*rcauchy(n, 0, 1) -
      ((B == 0)*qnorm(tau, 0, 1) + (B == 1)*qcauchy(tau, 0, 1))  ## mixture 0.9*N(0,1) + 0.1 Cauchy(0,1)
  }
  if (modtype == 1){## iid  case
    err <- e
  }  else if(modtype==2){ ## heteroscedasticity case
    err <- e*(1 - 0.2*x)  
  }
  y <- as.vector(xz %*% beta00) + err
  ##return(data.frame(y=y, x=x, z=z))
  return(cbind(y, x, z, rep(1,n)))
}


## CUSUM statistic
############################################
## tau: quantile
CUSUM.test <- function(x,y,z,tau){
  n <- length(y)
  datnew <- data.frame( y=y, x=x, z=z )
  fit <- rq(y~x+z, tau=tau, data=datnew, method="br")
  bet <- fit$coefficients
  res <- fit$residuals
  
  deri.checkfun <- function(u){
    tau-ifelse(u<0, 1, 0)
  }
  
  tt <- seq(min(x) + 0.01, max(x) - 0.01, length = 100)
  Rn <- c(NA, length(tt))
  for(k in 1:length(tt)){
    Rn[k] <- 1/sqrt(n) *sum(deri.checkfun(res)*(x-tt[k])*ifelse(x<tt[k], 1, 0))
  }
  L <- which.max(abs(Rn))
  return(c(max(abs(Rn)), tt[L]))
}

#CUSUM.test(x,y,z,tau)


############################################
## n: sample size
## M: the number of wild interval
random.intervals <-	function(n,M){
  n <- as.integer(n)
  M <- as.integer(M)
  intervals <- matrix(0,nrow=M,ncol=2)
  intervals[,1] <- ceiling(runif(M)*(n-1))
  intervals[,2] <- intervals[,1]+ ceiling(runif(M)*(n-intervals[,1]))
  intervals
}

#Inter <- random.intervals(n,M)

############################################
## M: the number of wild interval
wbs.quantile <- function(x, y, z, tau, M){
  data <- cbind(y,x,z) 
  n <- length(x)
  Inter <- random.intervals(n,M)
  #thre.para <- C*sqrt(2*n)
  
  cusum.break <- matrix(NA, M, 2)   ## save c(cusum, threshold)
  for(i in 1:M){
    dat1 <- data[Inter[i,1]: Inter[i,2],]
    y1 <- dat1[,1]
    x1 <- dat1[,2]
    z1 <- dat1[,3]
    tt <- seq(min(x1) + 0.01, max(x1) - 0.01, length = 100)
    #cusum.break[i,] <- CUSUM.test(x1, y1, z1, tau)
    test = try(CUSUM.test(x1, y1, z1, tau), silent=TRUE)
    if ('try-error' %in% class(test)){
      cusum.break[i,] <- c(0, NA)
    }else {
      cusum.break[i,] <- test  #CUSUM.test(x1, y1, z1, tau)
    }
  }
  local=which.max(cusum.break[,1])
  return(cusum.break[local, ])
}


#data1 <- cusum.break[order(cusum.break[,2]),] 
#wbs.quantile(x,y,z,tau,M)
#cp.break <- wbs.quantile(x,y,z,tau,M)[2]

############################################
## penalty function
## cp.break: the change points 
## kapa: parameter in SSIC function
penalty.ssci <- function(data, cp.break, kapa){
  n <- nrow(data)
  x <- data[,2]
  y <- data[,1]
  z <- data[,3]
  CP <- sort(cp.break)
  num <- length(cp.break)
  if(num==0){
    p <- 2+ncol(as.matrix(z,n,length(z)/n))
    datnew <- data.frame( y=y, x=x, z=z)
    fit <- rq(y~x+z, tau=tau, data=datnew, method="br")
    hat.y <- y-cbind(data[,4],data[,2],data[,3])%*%fit$coefficients
    loss <- sum(hat.y*(tau-ifelse(hat.y<0, 1, 0)))/n
    SIC <- log(loss)+log(n)/2/n*p
    sSIC <- n/2*log(loss)+num*log(n)^kapa
  } else{
    XX <- matrix(NA, n, num)
    p <- 2+ncol(as.matrix(z,n,length(z)/n))+2*num
    for(i in 1:num){
      XX[,i] <- pmax(x-CP[i],0)
    }
    datnew <- data.frame( y=y, x=x, XX, z=z )
    fit <- rq(y~x+XX+z, tau=tau, data=datnew, method="br")
    hat.y <- y-cbind(data[,4],data[,2],XX, data[,3])%*%fit$coefficients
    loss <- sum(hat.y*(tau-ifelse(hat.y<0, 1, 0)))/n
    SIC <- log(loss)+log(n)/2/n*p
    sSIC <- n/2*log(loss)+num*log(n)^kapa
  }
  return(sSIC)
}


#penalty.ssci(data, cp.break, kapa)
## wbs test procedure for quantile regression
##################################################
## Kmax: max number of change point
## 


CPtest.wbs.quan <- function (data, Kmax){
  result <- list()
  data <- data[order(data[,2]),] 
  x <- data[,2]
  y <- data[,1]
  z <- data[,3]
  CP.break <- CPs.cusum <- NULL
  int <- 0
  
  while(int < Kmax){
    int <- int+1
    
    if(length(CP.break)==0){
      CP.new1 <- wbs.quantile(x, y, z, tau, M)
      CPs.cusum <- CP.new1
      CP.new <- CP.new1[2]
    } else{
      CP.break <- as.vector(sort(CP.break))
      CP.cusum <- matrix(0, length(CP.break)+1, 2)
      
      for(i in 1:(length(CP.break)+1)){
        
        if(i==1){
          index <- which(x <= CP.break[i])
        } else if(i>1 & i<=length(CP.break)){
          index <- which(x <= CP.break[i] & x>CP.break[i-1])
        } else if(i==length(CP.break)+1){
          index <- which( x > CP.break[i-1])
        }
        
        xx <- x[index]
        yy <- y[index]
        zz <- z[index]
        #print(head(xx))
        CP.cusum[i,] <-  wbs.quantile(xx,yy,zz,tau,M) 
      }
      
      add.id <- which.max(CP.cusum[,1])
      CP.add <- CP.cusum[add.id, ]
      CP.new <- c( CP.break, CP.add[2])
      CPs.cusum <- rbind(CPs.cusum,  CP.add)
      if(min(abs(CP.break-CP.add[2]))<0.5) break
    }
    
    a <- try(penalty.ssci(data, CP.new, kapa), silent=TRUE)
    if ('try-error' %in% class(a)){
      int <- int-1
      pena.new <- -Inf
    } else{
      pena.new <- a
    }
    pena.break <- penalty.ssci(data, CP.break, kapa)
    if(pena.break <= pena.new) break
    
    CP.break <-  CP.new
  }    
  
  result$CP.break <- CP.break
  result$CPs.cusum <- CPs.cusum
  return(result)
}

#CPtest.wbs.quan(data, Kmax)


################################################################
## main run 
est.number<- rep(NA, NS)

for(i in 1:NS){
  dat <- bentdat(n, beta00, t00, tau, modtype, etype)
  data <- dat[order(dat[,2]),] 
  test.result <- CPtest.wbs.quan(data, Kmax)
  if ('try-error' %in% class(test.result)){
    est.number[i] <- NA
  } else{
    est.number[i] <- length(test.result$CP.break) 
  }
}


  





