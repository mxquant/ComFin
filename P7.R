# Explicit Finite-Difference Method 
efd <- function(s0, k, r, sig, t, dt, ds, log, put, amer){
  
  if (log) {  # if it is change in log price (dX)
    n <- floor((log(s0)-log(s0-3*sig*s0))/ds)
    # defind M and N
    M <- t/dt+1
    N <- 2*n+1

    dx <- ds
    # terminal sT
    sT <- exp(c(log(s0)+dx*(n:-n)))
  
    pu <- dt*(sig^2/(2*dx^2) + (r-sig^2/2)/(2*dx))
    pm <- 1 - dt*sig^2/dx^2 - r*dt
    pd <- dt*(sig^2/(2*dx^2) - (r-sig^2/2)/(2*dx))
    # A matrix
    A <- diag(pm, N-2, N-2)
    diag(A[-1,]) <- pu
    diag(A[,-1]) <- pd
    # set offset constans
    b <- c(pu, pd)
  } else { # if it is (dS)
    # defind M and N
    M <- t/dt+1
    sT <- seq(2*s0, 0, -ds)
    N <- length(sT)
    
    # calculate pu, pm, pd
    j <- (N-2):1
    pu <- 0.5*dt*(sig^2*j^2 + r*j)
    pm <- 1-dt*(sig^2*j^2+r)
    pd <- 0.5*dt*(sig^2*j^2 - r*j)
    # set A
    A <- diag(pm)
    diag(A[-1,]) <- pu[2:length(pu)]
    diag(A[,-1]) <- pd[1:length(pd)-1]
    # set offset constans
    b <- c(pu[1], pd[length(pd)])
  }

  payoff <- diag(0, N, M)
  if (put){
    payoff[,M] <- pmax(k-sT,0)
    payoff[1,] <- 0
    payoff[N,] <- (k-min(sT))*exp(-r*seq(t,0,-dt))
  } 
  else {
    payoff[,M] <- pmax(sT-k,0)
    payoff[1,] <- (max(sT)-k)*exp(-r*seq(t,0,-dt))
    payoff[N,] <- 0
  } 

  #loop through M
  for (i in (M-1):1){
    if (amer){  # if American option
      cv <- payoff[,i]
      cv[2:(N-1)] <- A%*%payoff[2:(N-1),i+1]
      cv[2] <- payoff[2,i]+b[1]*payoff[1,i+1]
      cv[N-1] <- payoff[N-1,i]+b[2]*payoff[N,i+1]
      payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
    } else { 
      payoff[2:(N-1),i] <- A%*%payoff[2:(N-1),i+1]
      payoff[2,i] <- payoff[2,i]+b[1]*payoff[1,i+1]
      payoff[N-1,i] <- payoff[N-1,i]+b[2]*payoff[N,i+1]
    }
  }
  return(payoff[floor(N/2)+1,1])
}

# Implicit Finite-Difference Method
ifd <- function(s0, k, r, sig, t, dt, ds, log, put, amer){
  if (log) {  #if dx
    n <- floor(log(s0)/ds)
    # defind M and N
    M <- t/dt+1
    N <- 2*n+1
    # define dx
    dx <- ds
    # terminal sT
    sT <- exp(c(log(s0)+dx*(n:-n)))

    pu = -0.5*dt*(sig^2/dx^2 + (r-sig^2/2)/dx)
    pm = 1 + dt*sig^2/dx^2 + r*dt
    pd = -0.5*dt*(sig^2/dx^2 - (r-sig^2/2)/dx)
    # A matrix
    A <- diag(pm, N-2, N-2)
    diag(A[-1,]) <- pu
    diag(A[,-1]) <- pd
    A <- rbind(c(1, -1, rep(0,N-2)), 
               cbind(c(pu, rep(0,N-3)),A,c(rep(0,N-3), pd)), 
               c(rep(0,N-2), 1, -1))
  } else { # if change in price
    # defind M and N
    M <- t/dt+1
    sT <- seq(2*s0, 0, -ds)
    N <- length(sT)
    
    j <- (N-2):1
    pu <- -0.5*dt*(sig^2*j^2 + r*j)
    pm <- 1+dt*(sig^2*j^2+r)
    pd <- 0.5*dt*(-sig^2*j^2 + r*j)
    # A matrix
    A <- diag(pm)
    diag(A[-1,]) <- pu[2:length(pu)]
    diag(A[,-1]) <- pd[1:length(pd)-1]
    # set offset constans
    b <- c(-pu[1], -pd[length(pd)])
  }

  payoff <- diag(0, N, M)
  if (put){
    payoff[,M] <- pmax(k-sT,0)
    payoff[1,] <- 0
    payoff[N,] <- (k-min(sT))*exp(-r*seq(t,0,-dt))
  } 
  else {
    payoff[,M] <- pmax(sT-k,0)
    payoff[1,] <- (max(sT)-k)*exp(-r*seq(t,0,-dt))
    payoff[N,] <- 0
  } 

  # find A inverse
  Ainv <- solve(A)
  # find option price
  if (log) {
    for (i in (M-1):1){
      # construct matrix B
      B <- payoff[,i+1]
      if (amer){  # American option
        cv <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      } else {  #european option
        payoff[,i] <- Ainv%*%B
      }
    }
  } else {
    # loop through M
    for (i in (M-1):1){
      # construct matrix B
      B <- payoff[2:(N-1),i+1]
      B[1] <- B[1] + payoff[1,i+1]*b[1]
      B[length(N)] <- B[length(N)] + payoff[N,i+1]*b[2]
      
      if (amer){ # if American option
        cv <- payoff[,i]
        cv[2:(N-1)] <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      } else { 
        payoff[2:(N-1),i] <- Ainv%*%B
      }
    }
  }
  return(payoff[floor(N/2)+1,1])
}


# Crank-Nicolson Finite-Difference Method
cnfd <- function(s0, k, r, sig, t, dt, ds, log, put, amer){
  if (log) { # dx
    n <- floor(log(s0)/ds)
    M <- t/dt+1
    N <- 2*n+1

    dx <- ds
    sT <- exp(c(log(s0)+dx*(n:-n)))

    pu = -0.25*dt*(sig^2/dx^2 + (r-sig^2/2)/dx)
    pm = 1 + dt*sig^2/(2*dx^2) + r*dt/2
    pd = -0.25*dt*(sig^2/dx^2 - (r-sig^2/2)/dx)
    # construct matrix A
    A <- diag(pm, N-2, N-2)
    diag(A[-1,]) <- pu
    diag(A[,-1]) <- pd
    A <- rbind(c(1, -1, rep(0,N-2)), 
               cbind(c(pu, rep(0,N-3)),A,c(rep(0,N-3), pd)), 
               c(rep(0,N-2), 1, -1))
  } else { # if ds

    M <- t/dt+1
    sT <- seq(2*s0, 0, -ds)
    N <- length(sT)

    j <- (N-2):1
    pu <- 0.25*dt*(sig^2*j^2 + r*j)
    pm <- -0.5*dt*(sig^2*j^2+r)
    pd <- 0.25*dt*(sig^2*j^2 - r*j)
    # set C
    A <- diag(1-pm)
    diag(A[-1,]) <- -pu[2:length(pu)]
    diag(A[,-1]) <- -pd[1:length(pd)-1]
    # set D
    D <- diag(1+pm)
    diag(D[-1,]) <- pu[2:length(pu)]
    diag(D[,-1]) <- pd[1:length(pd)-1]
    # set offset constans
    b <- c(pu[1], pd[length(pd)])
  }

  payoff <- diag(0, N, M)
  if (!put){ #call
    payoff[,M] <- pmax(sT-k,0)
    payoff[1,] <- (max(sT)-k)*exp(-r*seq(t,0,-dt))
    payoff[N,] <- 0
  } 
  else { #put
    payoff[,M] <- pmax(k-sT, 0)
    payoff[1,] <- 0
    payoff[N,] <- (k-min(sT))*exp(-r*seq(t,0,-dt))
  } 

  # find A inverse
  Ainv <- solve(A)
  # find option price
  if (log) {
    for (i in (M-1):1){
      # construct matrix B
      B <- payoff[,i+1]
      B[2:(N-1)] <- -pu*payoff[1:(N-2),i+1] - (pm-2)*payoff[2:(N-1),i+1] - pd*payoff[3:N,i+1] 
      # if Europian option
      if (!amer){
        payoff[,i] <- Ainv%*%B
      } else { 
        cv <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      }
    }
  } else {
    # loop through M
    for (i in (M-1):1){
      # construct matrix B
      B <- D%*%payoff[2:(N-1),i+1]
      B[1] <- B[1] + payoff[1,i+1]*b[1] + payoff[1,i]*b[1]
      B[length(N)] <- B[length(N)] + payoff[N,i+1]*b[2] + payoff[N,i+1]*b[2]
      # if Europian option
      if (!amer){
        payoff[2:(N-1),i] <- Ainv%*%B
      } else { # if American option
        cv <- payoff[,i]
        cv[2:(N-1)] <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      }
    }
  }
  return(payoff[floor(N/2)+1,1])
}

# Black-schole pricing model for European options
bsCall <- function(s0,t,x,r,sigma){
  d1 <- (log(s0/x)+(r+sigma^2/2)*t)/(sigma*sqrt(t))
  d2 <- d1-sigma*sqrt(t)
  c <- s0*pnorm(d1)-x*exp(-r*t)*pnorm(d2) 
  return(c)
}
bsPut <- function(s0,t,x,r,sigma){
  d1 <- (log(s0/x)+(r+sigma^2/2)*t)/(sigma*sqrt(t))
  d2 <- d1-sigma*sqrt(t)
  p <- x*exp(-r*t)*pnorm(-d2)-s0*pnorm(-d1) 
  return(p)
}


### Q1 ###
s0 <- 10
k <- 10
r <- 0.04
sig <- 0.2
t <- 0.5
dt <- 0.002
dx <- c(sig*sqrt(dt),sig*sqrt(3*dt),sig*sqrt(4*dt))

# i. European Put
p1=p2=p3=numeric(3)
# (a) Explicit Finite-Difference Method
p1 = sapply(dx,FUN=efd,s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=T,put=T,amer=F)
# (b) Implicit Finite-Difference Method
p2 = sapply(dx,FUN=ifd,s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=T,put=T,amer=F)
# (c) Crank-Nicolson Finite-Difference Method
p3 <- sapply(dx,FUN=cnfd,s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=T,put=T,amer=F)
# output
cbind(EFD = c(dt=p1[1],"3dt"=p1[2],"4dt"=p1[3]), IFD = p2, CNFE = p3)

# ii. Comparison
s0 = 4:16
pbs = sapply(s0, FUN=bsPut, t=t, x=k, r=r, sigma=sig)
for (i in 1:3){
  efd_put <- sapply(s0,FUN=efd, k=k, r=r, sig=sig, t=t, dt=dt, ds=dx[i], log=T,put=T,amer=F)
  ifd_put <- sapply(s0,FUN=ifd, k=k, r=r, sig=sig, t=t, dt=dt, ds=dx[i], log=T,put=T,amer=F)
  cnfd_put <- sapply(s0,FUN=cnfd, k=k, r=r, sig=sig, t=t, dt=dt, ds=dx[i], log=T,put=T,amer=F)

  err <- cbind(efd_put, ifd_put, cnfd_put)-pbs

  matplot(s0, abs(err), type ="l", lwd = 2, lty=1, main ="Error against Black-Scholes")
  legend("topleft", legend = colnames(err), lwd=1, col = 1:3)
}


### Q2 ###
s0 <- 10
k <- 10
r <- 0.04
sig <- 0.2
t <- 0.5
dt <- 0.002
ds <- c(0.5,1,1.5)

# i. American
# (a) Explicit Finite-Difference 
Ca <- sapply(ds,FUN=efd, s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=F, put=F, amer=T)
Pa <- sapply(ds,FUN=efd, s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=F, put=T, amer=T)
# (b) Implicit Finite-Difference 
Cb <- sapply(ds,FUN=ifd, s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=F, put=F, amer=T)
Pb <- sapply(ds,FUN=ifd, s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=F, put=T, amer=T)
# (c) Crank-Nicolson 
Cc <- sapply(ds,FUN=cnfd, s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=F, put=F, amer=T)
Pc <- sapply(ds,FUN=cnfd, s0=s0, k=k, r=r, sig=sig, t=t, dt=dt, log=F, put=T, amer=T)
# output
cbind(EFD = c("0.5ds"=Ca[1],"1.0ds"=Ca[2],"1.5ds"=Ca[3]), IFD = Cb, CNFE = Cc)  #calls
cbind(EFD = c("0.5ds"=Pa[1],"1.0ds"=Pa[2],"1.5ds"=Pa[3]), IFD = Pb, CNFE = Pc)  #puts

# ii. Graphs
s0 <- 4:16
cmat=rep(0,length(s0))
pmat=rep(0,length(s0))
for (i in 1:3){
  # calls
  efd_c <- sapply(s0,FUN=efd, k=k, r=r, sig=sig, t=t, dt=dt, ds=ds[i],log=F, put=F, amer=T)
  ifd_c <- sapply(s0,FUN=ifd, k=k, r=r, sig=sig, t=t, dt=dt, ds=ds[i],log=F, put=F, amer=T)
  cnfd_c <- sapply(s0,FUN=cnfd, k=k, r=r, sig=sig, t=t, dt=dt, ds=ds[i],log=F, put=F, amer=T)
  # puts
  efd_p <- sapply(s0,FUN=efd, k=k, r=r, sig=sig, t=t, dt=dt, ds=ds[i],log=F, put=T, amer=T)
  ifd_p <- sapply(s0,FUN=ifd, k=k, r=r, sig=sig, t=t, dt=dt, ds=ds[i],log=F, put=T, amer=T)
  cnfd_p <- sapply(s0,FUN=cnfd, k=k, r=r, sig=sig, t=t, dt=dt, ds=ds[i],log=F, put=T, amer=T)
  
  cmat=cbind(cmat,efd_c, ifd_c, cnfd_c)
  pmat=cbind(pmat,efd_p, ifd_p, cnfd_p)
}
matplot(s0, cmat[,-1], type ="b", pch=(1:9),lwd = 1.5, main="Call Prices",xlab="Stock price",ylab="Call Price")
legend("topleft", legend = c("EDF 0.5ds", "IDF 0.5ds", "CNDF 0.5ds","EDF 1.0ds", "IDF 1.0ds", "CNDF 1.0ds",
                             "EDF 1.5ds", "IDF 1.5ds", "CNDF 1.5ds"), pch=(1:9),lwd = 1.5,col = 1:9)

matplot(s0, pmat[,-1], type ="b", pch=(1:9),lwd = 1.5,main="Put Prices ", xlab = "Stock price", ylab = "Put Price")
legend("topright", legend = c("EDF 0.5ds", "IDF 0.5ds", "CNDF 0.5ds","EDF 1.0ds", "IDF 1.0ds", "CNDF 1.0ds",
                              "EDF 1.5ds", "IDF 1.5ds", "CNDF 1.5ds"), pch=(1:9),lwd = 1.5,col = 1:9)
