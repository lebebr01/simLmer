##Function that simulates errors
err.sim <- function(errorVar, n, p, serCor, serCorVal, err.dist){
  require(VGAM)
  require(MASS)
  if(serCor == "ARMA" & length(serCorVal) < 2) stop("Incorrect dimensions serCorVal")
  if(err.dist == "norm"){
    if(serCor == "AR"){
      err <- unlist(lapply(1:n,function(x){arima.sim(serCorVal,p, sd = sqrt(errorVar))}))
    } else {
      if(serCor == "MA"){
        err <- unlist(lapply(1:n,function(x){arima.sim(serCorVal,p, sd = sqrt(errorVar))}))
      } else {
        if(serCor == "ARMA"){
          err <- unlist(lapply(1:n,function(x){arima.sim(serCorVal,p, sd = sqrt(errorVar))}))
        } else {
          # generate multivariate normal error terms with zero mean 
          d2 <- (errorVar)*diag(p) 
          err <- matrix(c(mvrnorm(n,rep(0,p),d2)) ,nrow=n*p, ncol = 1)
        }
      }
    }
  } else {
    if(err.dist == "lap"){
      err <- unlist(lapply(1:n, function(x){rlaplace(p,0,1)*chol((errorVar/2))}))
    } else {
      if(err.dist == "chi"){
        err <- unlist(lapply(1:n, function(x){((rchisq(p,1)-1)/sqrt(2))*sqrt(errorVar)}))
      } else {
        #Note this does bimodal distribution with mean 0 and variance approx .64.
        #Does not change with errorVar
        err <- unlist(lapply(1:n, function(x){ c(rnorm(p/2,mean=.7,sd=.4),rnorm(p/2,mean=-.7,sd=.4))}))
      }
    }
  }
  err
}