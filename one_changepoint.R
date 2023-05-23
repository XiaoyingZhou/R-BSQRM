#############################################################################################
###   add kernel smooth
###   three methods: gria search, estimation by Das(2016) method, estimation by kernel method
###   f(0) estimated by the difference quotient
###   
#############################################################################################

#### one change point
#rm(list=ls())
ptm <- proc.time()
#set.seed(2017)
library(quantreg)

####  simulated data                                               
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
    e <- rcauchy(n, 0, 0.1) - qcauchy(tau, 0, 0.1)   ## Cauchy distribution
  } else if(etype ==5){
    prop <- 0.1
    B <- rbinom(n, 1, prop)
    e <-(B == 0)*rnorm(n, 0, 0.5) + (B == 1)*rcauchy(n, 0, 0.1) -
      ((B == 0)*qnorm(tau, 0, 0.5) + (B == 1)*qcauchy(tau, 0, 0.1))  ## mixture 0.9*N(0,1) + 0.1 Cauchy(0,1)
  }
  if (modtype == 1){## iid  case
    err <- e
  }  else if(modtype==2){ ## heteroscedasticity case
    err <- e*(1 + 0.2*x)  
  }
  y <- as.vector(xz %*% beta00) + err
  ##return(data.frame(y=y, x=x, z=z))
  return(cbind(y, x, z, rep(1,n)))
}


###################################################################
##  grid search algorithm                                        ##
##                                                               ##
###################################################################
## check function

checkfun <- function(u){
  u*(tau-ifelse(u<0, 1, 0))
}

qnregprof <- function(y, x, z) 
{
  tt <- seq(quantile(x,0.1), quantile(x,0.9), length=100)   # grid search interval
  sst<- NULL
  bet <- matrix(0, length(tt), 4)
  for (kk in 1:length(tt)){
    datnew<- data.frame(
      y=y, x1=x, x2=pmax(x-tt[kk],0), z= z
    )
    fit<- rq(y~x1+x2+z, tau=tau, data= datnew, method="br")
    bet[kk, ] <- fit$coef
    sst[kk] <- sum(checkfun(fit$residuals))
  }
  bp.grid <- min(tt[sst== min(sst)])
  bet.grid <- bet[tt==bp.grid, ]
  return(list(bet= bet.grid, bp=bp.grid))
}

Ku = function(u)
{
  # Epanechnikov kernel
  3/4*(abs(u)<=1)*(1-u^2)
}

qnregEse <- function(y, x, z, bet.est, t.est) 
{
  n <- length(y)
  alpha <- bet.est[1]
  beta1 <- bet.est[2]
  beta2 <- bet.est[3]
  gamm <- bet.est[4]
  gest <- alpha + beta1 * x + beta2 *pmax(x-t.est, 0) + gamm * z 
  res <- y - gest  
  hn <- 1.06* n^(-1/5) * sd(res)  ## the bandwidth by Silverman's rule of thumb
  ker <- Ku(res/hn)/hn
  
  gdev <- cbind(1, x, pmax(x-t.est, 0), z, (-1)*beta2*ifelse(x>t.est, 1, 0))
  Jn <- tau*(1-tau)/n * t(gdev)%*% gdev 
  Ln <- n^(-1) * t(gdev) %*% diag(ker) %*% gdev 
  
  Ln.inv <- solve(Ln+1e-8)
  Sig_thet <- n^(-1) * Ln.inv %*% Jn %*% Ln.inv
  ese <- sqrt(diag(Sig_thet))
  return(ese)
}



###############################################################################
##  estimation by Das(2016) method, f(0) estimated by the difference quotient
##                                                               
###############################################################################

fit.func <- function(theta){
  q <- rep(NA, n) ## save estimate value of q function 
  der.q <- rep(NA, n) ## save estimate value of derive of q function
  for(i in 1:n){
    if(x[i]<theta[5]-hn){
      q[i] <- 0
      der.q[i] <- 0
    }else if(theta[5]-hn<=x[i]&&x[i]<=theta[5]+hn){
      q[i] <- (x[i]-theta[5]+hn)^2/4/hn
      der.q[i] <- -2*(x[i]-theta[5]+hn)/4/hn
    }else {
      q[i] <- x[i]-theta[5]
      der.q[i] <- -1 
    }
  }
  return(rbind(q,der.q))
}


est.das <- function(hn, tau, theta0){
  H <- 0
  loss.fun <- function(tau,theta){
    for(i in 1:n){
      H1 <- 0
      if(x[i]<theta[5]-hn){
        H1 <- H1+theta[3]*0
      }else if(theta[5]-hn<=x[i]&&x[i]<=theta[5]+hn){
        H1 <- H1+theta[3]*(x[i]-theta[5]+hn)^2/4/hn
      }else {
        H1 <- H1+theta[3]*(x[i]-theta[5])
      }

      if(y[i]<theta[1]+theta[2]*x[i]+H1+theta[4]*z[i]){
        H <- H+(tau-1)*(y[i]-theta[1]-theta[2]*x[i]-theta[4]*z[i]-H1)
      }else {
        H <- H+tau*(y[i]-theta[1]-theta[2]*x[i]-theta[4]*z[i]-H1)
      }  
    }
    return(H)
  }
  
  fit <- optim(par=theta0, fn=function(theta){loss.fun(tau,theta)}, method = "BFGS")
  theta11 <- fit$par
  fit1 <- optim(par=theta0, fn=function(theta){loss.fun(tau-hN,theta)}, method = "BFGS")
  theta1 <- fit1$par
  fit2 <- optim(par=theta0, fn=function(theta){loss.fun(tau+hN,theta)}, method = "BFGS")
  theta2 <- fit2$par
  
  Q <- fit.func(theta11)
  q <- Q[1,]
  der.q <- Q[2,]
  
  Q1 <- fit.func(theta1)
  q1 <- Q1[1,]
  Dat1 <- cbind(1,x,q1, z) 
  f1 <- Dat1%*% theta1[-5]
  
  Q2 <- fit.func(theta2)
  q2 <- Q2[1,]
  Dat2 <-  cbind(1,x,q2, z)  ## cbind(data[,2], data[,3], q2, data[,5])
  f2 <- Dat2%*% theta2[-5]
  
  f.error <- 2*hN/(f2-f1)
  C <- 0
  V <- 0
  for(i in 1:n){
    h <- c(1,x[i], q[i], z[i], theta11[3]*der.q[i]) ## 
    C <- C+h%*%t(h)/n
    V <- V-f.error[i]*h%*%t(h)/n
  }
  Var <- tau*(1-tau)*solve(V)%*%C%*%t(solve(V))
  return(rbind(fit$par,sqrt(diag(Var)/n)))  
}



###############################################################################
##  estimation by kernel method, f(0) estimated by the difference quotient
##                                                               
###############################################################################

fit.func.ker <- function(theta){
  q <- (x-theta[5])*pnorm((x-theta[5])/hn.ker) ##  estimate value of kernel function 
  return(q)
}


est.Kern <- function(hn.ker, tau, theta0){
  H <- 0
  loss.fun = function(tau,theta){
    for (i in 1:n) {
      if(y[i] < theta[1]+theta[2]*x[i]+theta[3]*(x[i]-theta[5])*pnorm((x[i]-theta[5])/hn.ker)+theta[4]*z[i]){
        H = H+(tau-1)*(y[i]-theta[1]-theta[2]*x[i]-theta[3]*(x[i]-theta[5])*pnorm((x[i]-theta[5])/hn.ker)-theta[4]*z[i])
      }
      else {H = H+tau*(y[i]-theta[1]-theta[2]*x[i]-theta[3]*(x[i]-theta[5])*pnorm((x[i]-theta[5])/hn.ker)-theta[4]*z[i])}
    }
    return(H)
  }
  
  fit <- optim(par=theta0, fn=function(theta){loss.fun(tau,theta)}, method = "BFGS")
  theta11 <- fit$par
  fit1 <- optim(par=theta0, fn=function(theta){loss.fun(tau-hN,theta)}, method = "BFGS")
  theta1 <- fit1$par
  fit2 <- optim(par=theta0, fn=function(theta){loss.fun(tau+hN,theta)}, method = "BFGS")
  theta2 <- fit2$par
  
  q <- fit.func.ker(theta11)
  
  q1 <- fit.func.ker(theta1)
  Dat1 <- cbind(1,x,q1, z) 
  f1 <- Dat1%*% theta1[-5]
  
  q2 <- fit.func.ker(theta2)
  Dat2 <-  cbind(1,x,q2, z)  ## cbind(data[,2], data[,3], q2, data[,5])
  f2 <- Dat2%*% theta2[-5]
  
  f.error <- 2*hN/(f2-f1)
  C <- 0
  V <- 0
  for(i in 1:n){
    # h:  estimate value of derive of kernel function
    h <- c(1, x[i], (x[i]-theta11[5])*pnorm((x[i]-theta11[5])/hn.ker), z[i], -theta11[3]*(pnorm((x[i]-theta11[5])/hn.ker)+(x[i]-theta11[5])*dnorm((x[i]-theta11[5])/hn.ker)/hn.ker)) 
    C <- C+h%*%t(h)/n
    V <- V-f.error[i]*h%*%t(h)/n
  }
  Var <- tau*(1-tau)*solve(V)%*%C%*%t(solve(V))
  return(rbind(fit$par,sqrt(diag(Var)/n)))  
}

#################################################################
### parameter seting

alp<- 0.05
z.alp<- qnorm(1-alp/2)
hn <- n^(-1)       ## bandwith in smooth function
hn.ker <- 0.2*n^(-1/5)  ## bandwith in Kernel smooth function
hN <- n^(-1/3)*z.alp^(2/3)*(1.5*(dnorm(qnorm(tau)))^2/(2*(qnorm(tau))^2+1))^(1/3) ## bandwith in kernel function 


## save for grid search
theta.gs <- matrix(NA, NS, 5)
ese.gs <- matrix(NA, NS, 5)
CI1.gs <- matrix(NA, NS, 5)
CI2.gs <- matrix(NA, NS, 5)

## save for Das(2016) method
theta.hat <- matrix(NA, NS, 5)
ese.hat <- matrix(NA, NS, 5)
CI1.hat <- matrix(NA, NS, 5)
CI2.hat <- matrix(NA, NS, 5)

## save for kernel method
theta.K <- matrix(NA, NS, 5)
ese.K <- matrix(NA, NS, 5)
CI1.K <- matrix(NA, NS, 5)
CI2.K <- matrix(NA, NS, 5)

for(i in 1:NS){
  data <- bentdat(n, beta00, t00, tau, modtype, etype)
  x <- data[,2]
  y <- data[,1]
  z <- data[,3]
  
  ## by grid search algorithm 
  fit.gs <- qnregprof(y, x, z)
  if(max(abs(c(fit.gs$bet, fit.gs$bp)-c(beta00,t00)))>2){
    theta.gs[i, ] <- c(NA,NA,NA,NA,NA)
    ese.gs[i, ] <- c(NA,NA,NA,NA,NA)
    CI1.gs[i, ] <- c(NA,NA,NA,NA,NA)
    CI2.gs[i, ] <- c(NA,NA,NA,NA,NA)
  }else {
    theta.gs[i, ] <- c(fit.gs$bet, fit.gs$bp )
    if(max(qnregEse(y, x, z, fit.gs$bet, fit.gs$bp))>1){
      ese.gs[i, ] <- c(NA,NA,NA,NA,NA)
    }else {
      ese.gs[i, ] <- qnregEse(y, x, z, fit.gs$bet, fit.gs$bp) 
    }
    ##ese.gs[i, ] <- qnregEse(y, x, z, fit.gs$bet, fit.gs$bp)
    CI1.gs[i, ] <- theta.gs[i, ] - z.alp*ese.gs[i, ]
    CI2.gs[i, ] <- theta.gs[i, ] + z.alp*ese.gs[i, ]
  }
  
  
  theta0 <- c(fit.gs$bet, fit.gs$bp)  ## initial value
  
  ## by Das(2016) method estimate
  est <- try(est.das(hn, tau, theta0), silent=TRUE)
  if('try-error' %in% class(est)){
    theta.hat[i,] <- c(NA,NA,NA,NA,NA)
    ese.hat[i,] <- c(NA,NA,NA,NA,NA)
    CI1.hat[i,] <- c(NA,NA,NA,NA,NA)
    CI2.hat[i,] <- c(NA,NA,NA,NA,NA)
  }else {
    if(max(abs(est[1,]-c(beta00,t00)))>2){
      theta.hat[i,] <- c(NA,NA,NA,NA,NA)
      ese.hat[i,] <- c(NA,NA,NA,NA,NA)
      CI1.hat[i,] <- c(NA,NA,NA,NA,NA)
      CI2.hat[i,] <- c(NA,NA,NA,NA,NA)
    }else {
      theta.hat[i,] <-  est[1,]
      if(max(est[2,])>1){
        ese.hat[i,] <- c(NA,NA,NA,NA,NA)
      }else {
        ese.hat[i,] <- est[2,] 
      }
      CI1.hat[i,] <- theta.hat[i,] - z.alp*ese.hat[i,]
      CI2.hat[i,] <- theta.hat[i,] + z.alp*ese.hat[i,]
    }
  }
  
  ## by kernel method
  est1 <- est.Kern(hn.ker, tau, theta0)
  if('try-error' %in% class(est1)){
    theta.K[i,]<-c(NA,NA,NA,NA,NA)
    ese.K[i,]<-c(NA,NA,NA,NA,NA)
    CI1.K[i,]<-c(NA,NA,NA,NA,NA)
    CI2.K[i,]<-c(NA,NA,NA,NA,NA)
  }else {
    if(max(abs(est1[1,]-c(beta00,t00)))>2){
      theta.K[i,]<-c(NA,NA,NA,NA,NA)
      ese.K[i,]<-c(NA,NA,NA,NA,NA)
      CI1.K[i,]<-c(NA,NA,NA,NA,NA)
      CI2.K[i,]<-c(NA,NA,NA,NA,NA)
    }else {
      theta.K[i,] <-  est1[1,]
      ese.K[i,] <- est1[2,]
      if(max(est1[2,])>1){
        ese.K[i,] <- c(NA,NA,NA,NA,NA)
      }else {
        ese.K[i,] <- est1[2,] 
      }
      CI1.K[i,] <- theta.K[i,] - z.alp*ese.K[i,]
      CI2.K[i,] <- theta.K[i,] + z.alp*ese.K[i,] 
    }
  }
}


true <- matrix(rep(c(beta00, t00),NS), nrow=NS, byrow=T)
## grid search
grid.est <- apply(theta.gs, 2, mean, na.rm=TRUE)
grid.bias <- grid.est-c(beta00,t00)
grid.sd <- apply(theta.gs, 2, sd, na.rm=TRUE)
grid.ese <- apply(ese.gs, 2, mean, na.rm=TRUE)
grid.MSE <- grid.bias^2+grid.sd^2
grid.CI1 <- apply(CI1.gs, 2, mean, na.rm=TRUE)
grid.CI2 <- apply(CI2.gs, 2, mean, na.rm=TRUE)
grid.AW <- grid.CI2 - grid.CI1
grid.CP <- apply( (CI1.gs< true)*(true<CI2.gs), 2, mean, na.rm=TRUE)

## Das(2016) method  estimate
pro.est <- apply(theta.hat, 2, mean, na.rm=TRUE)
pro.bias <- pro.est-c(beta00,t00)
pro.sd <- apply(theta.hat, 2, sd, na.rm=TRUE)
pro.ese <- apply(ese.hat, 2,  mean, na.rm=TRUE)
pro.MSE <- pro.bias^2+pro.sd^2
pro.CI1 <- apply(CI1.hat, 2,  mean, na.rm=TRUE)
pro.CI2 <- apply(CI2.hat, 2,  mean, na.rm=TRUE)
pro.AW <- pro.CI2 - pro.CI1
pro.CP <- apply( (CI1.hat< true)*(true<CI2.hat), 2, mean, na.rm=TRUE)


## kernel method estimate
kern.est <- apply(theta.K, 2, mean, na.rm=TRUE)
kern.bias <- kern.est-c(beta00,t00)
kern.sd <- apply(theta.K, 2, sd, na.rm=TRUE)
kern.ese <- apply(ese.K, 2,  mean, na.rm=TRUE)
kern.MSE <- kern.bias^2+kern.sd^2
kern.CI1 <- apply(CI1.K, 2,  mean, na.rm=TRUE)
kern.CI2 <- apply(CI2.K, 2,  mean, na.rm=TRUE)
kern.AW <- kern.CI2 - kern.CI1
kern.CP <- apply( (CI1.K< true)*(true<CI2.K), 2, mean, na.rm=TRUE)

out<-rbind(
  rbind(grid.est, grid.bias, grid.sd, grid.ese, grid.MSE, grid.CP, grid.AW, grid.CI1, grid.CI2), 
  rbind(pro.est, pro.bias, pro.sd, pro.ese, pro.MSE, pro.CP, pro.AW, pro.CI1, pro.CI2),
  rbind(kern.est, kern.bias, kern.sd, kern.ese, kern.MSE, kern.CP, kern.AW, kern.CI1, kern.CI2)
) 

colnames(out) = c("beta0", "beta1", "beta2", "gamma", "t")

folder="E:\\paper_simulation\\Das_BentReg\\bent line quantile\\AllNewRsults\\est1-t-0\\"
filename = paste("est1","etype", etype, "modtype", modtype,"tau", tau,"n", n, sep="-") 
write.csv(out, file=paste(folder, filename, ".csv", sep=''),  row.names=T)

proc.time()-ptm

