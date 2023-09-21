library(goffda)


simulate_LM <- function(lambda,sigma,alpha,lambda.t){
  
  current.optimum <- rnorm(1,mean=0,sd=sigma/sqrt(2*alpha))
  
  #initialise the state of the lineage
  pop1 <- list(phenotype=c(0),
              time=c(0),
              current.optimum=current.optimum,
              upper.neighbour=current.optimum+rexp(1,rate=lambda),
              lower.neighbour=current.optimum-rexp(1,rate=lambda))
           
  #10,000 years in increments of 100 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- goffda::r_ou(1,t=seq(0,10^-2,len=101),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  
  time <- pop1$time[length(pop1$time)]
  
  while(time < 111){
    idx.lower <- rexp(1,lambda.t)
    t0 <- pop1$time[length(pop1$time)]
    phenotype.t <- goffda::r_ou(1,t=seq(0,idx.lower,len=101),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
    pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
    pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
    
    if(runif(1) < 0.5){
      pop1$current.optimum <- pop1$lower.neighbour
      pop1$lower.neighbour <- pop1$current.optimum - rexp(1,rate=lambda)
      pop1$upper.neighbour <- pop1$current.optimum + rexp(1,rate=lambda)
    } else{
      pop1$current.optimum <- pop1$upper.neighbour
      pop1$lower.neighbour <- pop1$current.optimum - rexp(1,rate=lambda)
      pop1$upper.neighbour <- pop1$current.optimum + rexp(1,rate=lambda)
    }
    time <- pop1$time[length(pop1$time)]
  }
  return(pop1)
}

#try <- simulate_LM(10,20,100000,0.2)
#plot(log10(try$time),try$phenotype,type='l')

#at which points in time should the process be sampled
#corresponding to time (in years): 10^2, 10^3, 10^3.5, 10^4, 10^4.5, 10^5, 10^5.3, 10^6, 10^6.7, 10^7.3, 10^7.7, 10^7.9

simulate_variance <- function(x=c(seed,lambda,sigma,alpha,lambda.t)){
  
  set.seed(x[1])
  n <- 100 #number of simulations of the process to run
  samplepoints <- 10^c(-4,-3,-2.5,-2,-1.5,-1,-0.7,0,0.7,1.3,1.7,1.9)
  
  t <- simulate_LM(lambda=x[2],sigma=x[3],alpha=x[4],lambda.t=x[5])
  
  idx <- numeric()
  for(i in 1:length(samplepoints)) idx[i] <- which.min(abs(samplepoints[i]-t$time))
  t <- t$phenotype[idx]
  
  for(i in 1:n){
    tt <- simulate_LM(lambda=x[2],sigma=x[3],alpha=x[4],lambda.t=x[5])
    idx <- numeric()
    for(i in 1:length(samplepoints)) idx[i] <- which.min(abs(samplepoints[i]-tt$time))
    tt <- tt$phenotype[idx]
    t <- cbind(t,tt)
  }
  #calculate variances from simulated divergence data
  vars <- numeric()
  for(i in 1:length(t[,1])){
    dat.t <- c(t[i,])
    vars[i] <- log(var(dat.t))
  }
  times <- c(2,3,3.5,4,4.5,5,5.3,6,6.7,7.3,7.7,7.9)
  names(vars) <- times
  
  return(vars)
}

#try <- simulate_variance(x=c(42,20,7,7000,0.5))
#system.time(simulate_variance(x=c(42,2,10,70000)))
#plot(names(try),try)





