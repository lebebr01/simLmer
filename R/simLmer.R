##############################################

#Attempt at new function

############################
#fixed = one sided formula for fixed effects to be included in the simulation, include 1 for intercept
#random = one sided formula for random effects to be included in the simulation, must be a subset of fixed
#fixed.param = fixed effect parameter values (i.e. beta weights)  Must be same length as fixed.
#random.param = variance of random effects, must be same length as random.
#w.var = number of within cluster variables - including intercept
#cov.param = list of mean and variance of fixed effects besides intercept/time and interactions, must be same order as fixed
#n = cluster sample size
#p = within cluster sample size
#errorVar = Error variance
#randCor = Correlation between random effects
#rand.dist = random effect distribution = either "lap", "chi", "norm", "bimod"
#err.dist = within cluster error distribution = either "lap", "chi", "norm", "bimod"
#serCor = Serial correlation present = either "AR", "MA", "ARMA", or "ID"
#serCorVal = Serial correlation values, a list of values to pass on to arima.sim
#data.str = type of data, either cross or long  #Currently not supported
#####################################

sim.lmer <- function(fixed, random, fixed.param, random.param, w.var, cov.param, n, p, errorVar, randCor, rand.dist, err.dist, serCor, serCorVal) {

  if(randCor > 1 | randCor < -1) stop("cor out of range")

  require(MASS)
  require(VGAM)
  require(Matrix)

  formula <- paste(fixed,random, sep="+(")[2]
  formula <- gsub("$",")", formula)
  formula <- gsub(" ", "", formula)
  
  fixed.vars <- attr(terms(fixed),"term.labels")    ##Extracting fixed effect term labels
  rand.vars <- attr(terms(random),"term.labels")   ##Extracting random effect term labels

     if(length(rand.vars)+1 != length(random.param)) stop("Random lengths not equal")
     if({length(fixed.vars)+1} != {length(fixed.param)}) stop("Fixed lengths not equal")

   rand.eff <- rand.eff.sim(random.param, randCor, n, rand.dist)

    
       n.vars <- length(fixed.vars)
       n.int <- length(grep(":",fixed.vars))
        Xmat <- matrix(nrow=n*p,ncol = n.vars+1-n.int)
         if(n.int == 0){
          colnames(Xmat) <- c("Int", fixed.vars)
         } else {
           int.loc <- grep(":", fixed.vars)
           colnames(Xmat) <- c("Int",fixed.vars[-int.loc])
         }

 if(w.var == 1){
  Xmat[,1] <- rep.int(1,times = n*p)
  } else {
   if(w.var == 2){
    Xmat[,1] <- rep.int(1, times = n*p)
    Xmat[,2] <- rep.int((1:p)-1,times = n)
   } else { 
     Xmat[,1] <- rep.int(1, times = n*p)
     Xmat[,2] <- rep.int((1:p)-1,times = n)
      for(j in 3:w.var){
       Xmat[,j] <- rep.int(rnorm(p,mean=cov.param[[j-2]][1],sd=cov.param[[j-2]][2]), times = n)
      }
   }
 }

 if(n.int == 0){
  if(w.var != n.vars+1){
     for(l in (w.var+1):(n.vars+1)){
       Xmat[,l] <- rep(rnorm(n, mean = cov.param[[l-2]][1], sd=cov.param[[l-2]][2]), each = p)
     }
  } else {
    next;
  }
 } else {
     num.no.int <- n.vars - n.int                  
       if(w.var != num.no.int+1){
          for(k in (w.var+1):(num.no.int+1)){
            Xmat[,k] <- rep(rnorm(n, mean = cov.param[[k-2]][1], sd=cov.param[[k-2]][2]), each = p)
          }
       } else {
         next;
       }
     Xmat <- model.matrix(fixed, data.frame(Xmat))
     
 }

 if(ncol(rand.eff) == 2){            ##Adding random effects to design matrix
  re <- rep(rand.eff,each=p)
  Zmat <- data.frame(cbind(re[1:(n*p)],re[(n*p+1):(n*p*2)]))
   colnames(Zmat) <- c("b0", "b1")
 } else {
  Zmat$b0 <- data.frame(rep(rand.eff, each = p))
 }

 #if(serCor == "AR" | serCor == "MA" | serCor == "ARMA" & is.list(serCorVal) == "FALSE") {stop("Incorrect dimensions serCorVal")}
  err <- err.sim(errorVar, n, p, serCor, serCorVal, err.dist)

 sim.data <- data.lmer(Xmat, fixed.param, n, p, rand.eff, err)
  
 Xmat <- data.frame(Xmat,Zmat,err,sim.data)
 Xmat$ID <- rep(1:n, each = p)
 Xmat
}