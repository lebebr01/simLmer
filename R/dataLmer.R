##Function that simulates correlated data
dataLmer <- function(Xmat, beta, n, p, rand.eff, err) {
  if(ncol(rand.eff) == 2){
    #ZmatT <- matrix(c(t(rand.eff)))  #Convert random effects into column matrix
    #XmatSin <- matrix(c(rep(1,p),(1:p)-1),nrow=p,ncol=2)
    #XmatR <- blockDiagonal(XmatSin, n*p, 2*n)
    #XmatTime <- matrix(c((1:p)-1),nrow=p,ncol=1)
    #XmatR <- blockDiagonal(XmatTime, n*p, n)  #Generate Block Diagonal Matrix
    Fbeta <-(Xmat %*% beta)  #Simulate average growth curve
    #Reff <- (XmatR %*% ZmatT)  #Simulate random effect portion - Not Working for N greater than 25
    b0 <- matrix(rep(rand.eff[,1], each = p),ncol=1)
    b1 <- rep(rand.eff[,2], each = p) * rep((1:p)-1, n)
    Reff <- b0 + b1
    sim.data <- Fbeta + Reff + err  #Adding everything together
    #sim.data <- Fbeta + b0 + b1 + err
    sim.data <- cbind(Fbeta, Reff, sim.data)  
    colnames(sim.data) <- c("Fbeta", "Reff", "sim.data")
    sim.data
  } else {
    Fbeta <- (Xmat %*% beta) 
    Reff <- Zmat 
    sim.data <- Fbeta + Reff + err
    sim.data <- cbind(Fbeta, Reff, sim.data)
    colnames(sim.data) <- c("Fbeta", "Reff", "sim.data")
    sim.data
  }
}