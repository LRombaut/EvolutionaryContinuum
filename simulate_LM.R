library(goffda)

simulate_LM <- function(lambda,sigma,alpha){
  
  current.optimum <- rnorm(1,mean=0,sd=sigma/sqrt(2*alpha))
  
  #initialise the state of the lineage
  pop1 <- list(phenotype=numeric(),
              time=numeric(),
              current.optimum=current.optimum,
              upper.neighbour=current.optimum+rexp(1,rate=lambda),
              lower.neighbour=current.optimum-rexp(1,rate=lambda))
 
  ##initialize the first million years of the population's evolution
  #the first 1 years...
  phenotype.t <- r_ou(1,t=seq(0,0.000001,len=11),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=0)
  pop1$phenotype <- phenotype.t$data
  pop1$time <- phenotype.t$argvals
  #10 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- r_ou(1,t=seq(0,0.00001,len=11),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  #100 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- r_ou(1,t=seq(0,0.0001,len=11),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  #1,000 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- r_ou(1,t=seq(0,0.001,len=11),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  #10,000 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- r_ou(1,t=seq(0,0.01,len=11),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  #100,000 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- r_ou(1,t=seq(0,0.1,len=11),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  #1,000,000 years
  t0 <- pop1$time[length(pop1$time)]
  phenotype.t <- r_ou(1,t=seq(0,1,len=101),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
  pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
  pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
  
  #if and when does the population exceed the thresholds to the neighbouring adaptive zones?
  idx.upper <- min(which(pop1$phenotype > pop1$current.optimum + 0.5*(pop1$upper.neighbour-pop1$current.optimum)))
  idx.lower <- min(which(pop1$phenotype < pop1$current.optimum - 0.5*(pop1$current.optimum-pop1$lower.neighbour)))
  
  time <- pop1$time[length(pop1$time)]
  
  while(time < 101){ #time measured in millions of years
    #the population has remained within the bounds of the adaptive zone
    if(is.infinite(idx.lower) & is.infinite(idx.upper)){
      t0 <- pop1$time[length(pop1$time)]
      phenotype.t <- r_ou(1,t=seq(0,1,len=101),mu=pop1$current.optimum,alpha=alpha,sigma=sigma,x0=pop1$phenotype[length(pop1$phenotype)])
      pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
      pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
    }
    if(idx.upper<idx.lower){
      t0 <- pop1$time[idx.upper]
      pop1$phenotype <- pop1$phenotype[1:idx.upper]
      pop1$time <- pop1$time[1:idx.upper]
      phenotype.t <- r_ou(1,t=seq(0,1,len=101),mu=pop1$upper.neighbour,alpha=alpha,sigma=sigma,x0=pop1$phenotype[idx.upper])
      pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
      pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
      if(!any(phenotype.t$data[-1] < pop1$phenotype[idx.upper])){
        pop1$current.optimum <- pop1$upper.neighbour
        pop1$upper.neighbour <- pop1$current.optimum+rexp(1,rate=lambda)
        pop1$lower.neighbour <- pop1$current.optimum-rexp(1,rate=lambda)
      } else{
        pop1$lower.neighbour <- pop1$current.optimum
        pop1$current.optimum <- pop1$upper.neighbour
        pop1$upper.neighbour <- pop1$upper.neighbour+rexp(1,rate=lambda)
      }
    }
    if(idx.lower<idx.upper){
      t0 <- pop1$time[idx.lower]
      pop1$phenotype <- pop1$phenotype[1:idx.lower]
      pop1$time <- pop1$time[1:idx.lower]
      phenotype.t <- r_ou(1,t=seq(0,1,len=101),mu=pop1$lower.neighbour,alpha=alpha,sigma=sigma,x0=pop1$phenotype[idx.lower])
      pop1$phenotype <- c(pop1$phenotype,phenotype.t$data[-1])
      pop1$time <- c(pop1$time,phenotype.t$argvals[-1]+t0)
      if(!any(phenotype.t$data[-1] > pop1$phenotype[idx.lower])){
        pop1$current.optimum <- pop1$lower.neighbour
        pop1$upper.neighbour <- pop1$current.optimum+rexp(1,rate=lambda)
        pop1$lower.neighbour <- pop1$current.optimum-rexp(1,rate=lambda)
      } else{
        pop1$upper.neighbour <- pop1$current.optimum
        pop1$current.optimum <- pop1$lower.neighbour
        pop1$lower.neighbour <- pop1$lower.neighbour-rexp(1,rate=lambda)
      }
    }
    idx.upper <- min(which(phenotype.t$data[-1] > pop1$current.optimum + 0.5*(pop1$upper.neighbour-pop1$current.optimum)))
    idx.upper <- length(pop1$phenotype) - length(phenotype.t$data[-1]) + idx.upper
    idx.lower <- min(which(phenotype.t$data[-1] < pop1$current.optimum - 0.5*(pop1$current.optimum-pop1$lower.neighbour)))
    idx.lower <- length(pop1$phenotype) - length(phenotype.t$data[-1]) + idx.lower
    
    time <- pop1$time[length(pop1$time)]
  }
  return(pop1)
}


simulate_variance <- function(x=c(seed,lambda,sigma,alpha)){
  
  set.seed(x[1])
  n <- 100 #number of simulations of the process to run
  samplepoints <- 10^c(-4,-3,-2.5,-2,-1.5,-1,-0.7,0,0.7,1.3,1.7,1.9)
  
  t <- simulate_LM(lambda=x[2],sigma=x[3],alpha=x[4])
  
  idx <- numeric()
  for(i in 1:length(samplepoints)) idx[i] <- which.min(abs(samplepoints[i]-t$time))
  t <- t$phenotype[idx]
  
  for(i in 1:n){
    tt <- simulate_LM(lambda=x[2],sigma=x[3],alpha=x[4])
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

#system.time(simulate_variance(x=c(42,2,10,70000)))





