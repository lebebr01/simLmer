##Function to simulate random effects
randEffSim <- function(random.param, cor, n, dist = c("lap","chi","norm")){
  
  require(MASS)
  require(VGAM)
  if(dist == "lap"){ 
    if(length(random.param) == 2) {
      b0 <- rlaplace(n,0,1) # individual intercepts' deviation from fixed intercept
      b1 <- rlaplace(n,0,1) # individual slopes' deviation from fixed slope
      reff<-cbind(b0,b1)
      cov <- cor*sqrt(random.param[1]*random.param[2])
      c <- matrix(c(random.param[1]/2,cov/2,cov/2,random.param[2]/2),byrow=T,ncol=2)
      reff1<-reff%*%chol(c)
      reff1
    } else {
      b0 <- rlaplace(n,0,1)
      reff <- b0*chol(random.param/2)
      reff
    }
  } else {
    if(dist == "chi"){ 
      if(length(random.param) == 2) {
        b0 <- rchisq(n,1) # individual intercepts' deviation from fixed intercept
        b1 <- rchisq(n,1) # individual slopes' deviation from fixed slope
        reff<-cbind(b0-1,b1-1)
        cov <- cor*sqrt(random.param[1]*random.param[2])
        c <- matrix(c(random.param[1]/2,cov/2,cov/2,random.param[2]/2),byrow=T,ncol=2)
        reff1<-reff%*%chol(c)
        reff1
      } else {
        b0 <- rchisq(n,1)
        reff1 <- (b0-1)*chol(random.param/2)
        reff1
      }
    } else {
      if(dist == "bimod"){
        if(length(random.param) == 2) {
          #Note only does specific simulation variances
          b0 <- c(rnorm(n/2,mean=.63,sd=.3),rnorm(n/2,mean=-.63,sd=.3))
          b1 <- c(rnorm(n/2,mean=.21,sd=.08),rnorm(n/2,mean=-.21,sd=.08))
          reff <- cbind(b0,b1)
          reff
        } else {
          #Note only does specific simulation variances
          b0 <- c(rnorm(p/2,mean=.63,sd=.3),rnorm(p/2,mean=-.63,sd=.3))
          b0
        }
      } else { 
        if(length(random.param) == 2){
          cov <- cor*sqrt(random.param[1]*random.param[2]) # correlation between random effects
          d <- matrix(c(random.param[1],cov,cov,random.param[2]),nrow=2,ncol=2) #var-cov matrix of random effects
          ind <- mvrnorm(n,c(0,0),d) # generate bivariate normal random effects
          colnames(ind) <- c("b0","b1")
          ind
          #print(ind)
          
        } else {
          d <- matrix(c(random.param), nrow=1,ncol=1)
          ind <- mvrnorm(n,0,d)
          colnames(ind) <- "b0"
          ind
          #print(ind)
        }
      }}}}