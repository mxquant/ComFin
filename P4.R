#ComFin Project 4
rm(list = ls())
library(derivmkts) 
library(quantmod)
library(data.table)


gbm=function(s0=100, sig=0.2, r=0.06, t=1, steps=365, n=1000, d=0){ 
  st <- matrix(NA, nrow = n, ncol = steps)
  for (i in 1:n) {
    st[i, ] = s0 * exp(cumsum((r-d-0.5*sig^2)*(t/steps) + (sig*(sqrt(t/steps))*rnorm(steps))))
  }
  return(st)
}


#Question 1
binom=function(u,d,pu,s0,x,r,t,n,amr=F,put=F){ 
  pd=1-pu
  k=n:0 
  st=s0*u^k*d^(n-k)   #last period stock prices
  if(put){c=pmax(x-st,0)}  else{c=pmax(st-x,0)}  #last period option values
  
  for(j in n:1){
    st=st/u   #for each time step, move prices 1 period back
    for(i in 1:j){
      c[i]=(pd*c[i+1]+pu*c[i])*exp(-r*t/n)  #continuing values
      if(amr){c[i]=ifelse(put,max(c[i],x-st[i]),max(c[i],st[i]-x))}  #check exercise values
    }
  }
  return(c[1])
}

q1a <- function(r, sigma, t, n){
  dt = t/n
  c = .5*(exp(-r*dt)+exp((r+sigma^2)*dt))
  d = c - sqrt(c^2 - 1)
  u = 1/d
  p = (exp(r*dt) - d)/(u - d)
  return(c(u, d, p))
}
q1b <- function(r, sigma, t, n){
  dt = t/n
  p = .5
  u = exp(r*dt)*(1+sqrt(exp(sigma^2*dt)-1))
  d = exp(r*dt)*(1-sqrt(exp(sigma^2*dt)-1))
  return(c(u, d, p))
}
q1c <- function(r, sigma, t, n){
  dt = t/n
  p = .5
  u = exp((r-sigma^2/2)*dt+sigma*sqrt(dt))
  d = exp((r-sigma^2/2)*dt-sigma*sqrt(dt))
  return(c(u, d, p))
}
q1d <- function(r, sigma, t, n){
  dt = t/n
  u = exp(sigma*sqrt(dt))
  d = exp(-sigma*sqrt(dt))
  p = .5+.5*(((r-sigma^2/2)*sqrt(dt))/sigma)
  return(c(u, d, p))
}
# Output
s0=32
k=30
r=0.05
sig=0.24
t=0.5
n=c(10,20,40,80,100,200,500)
sol=matrix(0,length(n),4)
for (i in 1:length(n)){
  udp1 <- q1a(r, sig, t, n[i])
  sol[i,1] <- binom(udp1[1], udp1[2], udp1[3], s0, k, r, t, n[i])

  udp2 <- q1b(r, sig, t, n[i])
  sol[i,2] <- binom(udp2[1], udp2[2], udp2[3], s0, k, r, t, n[i])

  udp3 <- q1c(r, sig, t, n[i])
  sol[i,3] <- binom(udp3[1], udp3[2], udp3[3], s0, k, r, t, n[i])

  udp4 <- q1d(r, sig, t, n[i])
  sol[i,4] <- binom(udp4[1], udp4[2], udp4[3], s0, k, r, t, n[i])
}

matplot(n,sol, type ="l", lwd=2,lty=1,main="Convergence Rates",xlab="steps", ylab ="Call Price", col = c(1:4))
abline(h=bscall(s0,k,sig,r,t,0), lty = 5)
legend("topright", legend = paste0("(",letters[1:4],")"), lwd = 2, col = c(1:4))


#Question 2
getSymbols(Symbols = 'AMZN', src = "yahoo", from='2015-04-24', to='2020-04-24')
prices=Ad(AMZN)
ret=(shift(prices,type="lead")-prices)/prices
s0 <- as.numeric(prices[nrow(prices)])
sig<-sd(ret,na.rm = T)*sqrt(252)
r <- 0.02
k <- as.numeric(round(s0*1.1/10)*10)
t <- as.numeric(as.Date("2021-01-15") - as.Date("2020-04-24"))/365 
n <- 1000

# use binomial method from 1(a)
udp = q1a(r, sig, t, n)
amzncall = binom(udp[1], udp[2], udp[3], s0, k, r, t, n)
amzncall

# (b) the market price of the call option is 229 on April 24, 2020	
mkt <- 229
to_solve=function(sigma){
  binom(q1a(r,sigma,t,n)[1],q1a(r,sigma,t,n)[2],q1a(r,sigma,t,n)[3],s0,k,r,t,n) - mkt
}
uniroot(to_solve,c(0.1,2))$root


#Question 3
s0 <- seq(20, 80, 2)
k <- 50
r <- 0.03
sig <- 0.2
t <- 0.3846
n <- 1000

# (i) use binomial method from 1(a)
delta = numeric(length(s0))
udp = q1a(r, sig, t, n)
oricall=sapply(s0,FUN=binom,u=udp[1],d=udp[2],pu=udp[3],k,r,t,n)
for (i in 1:length(s0)){
  newcall = binom(udp[1], udp[2], udp[3], s0[i]+0.01, k, r, t, n)
  delta[i] <- (newcall - oricall[i])/0.01
}
plot(s0,delta, type = "l", lwd = 2,main = "Delta vs S0", xlab = "S0", ylab = "Delta")


# (iii)
theta <- numeric(length(s0))
newudp = q1a(r,sig,t-0.0001,n)
for (i in 1:length(s0)){
  newcall = binom(newudp[1], newudp[2], newudp[3], s0[i], k, r, t-0.0001, n)
  theta[i] = (newcall - oricall[i])/0.0001
}
plot(s0,theta, type = "l", lwd = 2, col = "red",main = "Theta vs S0", xlab = "S0", ylab = "Theta")

# (iv)
gamma = numeric(length(s0))
for (i in 1:length(s0)){
  c1 = binom(udp[1], udp[2], udp[3], s0[i]-1, k, r, t, n)
  c2 = binom(udp[1], udp[2], udp[3], s0[i]+1, k, r, t, n)
  gamma[i] <- (c1+c2-2*oricall[i])/(1^2)
}
plot(s0,gamma, type = "l", col="blue",lwd = 2,main = "Gamma vs S0", xlab = "S0", ylab = "Gamma")

# (v)
vega <- numeric(length(s0))
newudp = q1a(r,sig+0.0001,t,n)
for (i in 1:length(s0)){
  newcall = binom(newudp[1], newudp[2], newudp[3], s0[i], k, r, t, n)
  vega[i] = (newcall - oricall[i])/0.0001/100
}
plot(s0,vega, type = "l", lwd = 2, col = "green",main = "Vega vs S0", xlab = "S0", ylab = "Vega")

# (vi)
rho <- numeric(length(s0))
newudp = q1a(r+0.0001,sig,t,n)
for (i in 1:length(s0)){
  newcall = binom(newudp[1], newudp[2], newudp[3], s0[i], k, r+0.0001, t, n)
  rho[i] = (newcall - oricall[i])/0.0001/100
}
plot(s0,rho, type = "l", lwd = 2, col = "purple",main = "Rho vs S0", xlab = "S0", ylab = "Rho")

# (ii)
tt = seq(0, 0.3846, 0.01)
s0 = 49
delta2=numeric(length(tt))
for (i in 1:length(tt)){
  newudp = q1a(r, sig, tt[i], n)
  oricall = binom(newudp[1], newudp[2], newudp[3], s0, k, r, tt[i], n)
  newcall = binom(newudp[1], newudp[2], newudp[3], s0+0.01, k, r, tt[i], n)
  delta2[i] = (newcall-oricall)/(0.01)
}
plot(tt,delta2, type = "l", lwd = 2, col = "orange", main = "Delta vs Time to Maturity", 
     xlab = "Time to Maturity", ylab = "Delta")


#Question 4
t=1
r=0.05
sig=0.3
k=100
s0=seq(80,120,4)
n=1000
#binom=function(u,d,pu,s0,x,r,t,n,amr=F,put=F)   q1a=function(r,sigma,t,n)
udp=q1a(r,sig,t,n)
sol4=matrix(0,length(s0),2)
sol4[,1]=sapply(s0,FUN=binom,u=udp[1],d=udp[2],pu=udp[3],x=k,r=r,t=t,n=n,put=T)
sol4[,2]=sapply(s0,FUN=binom,u=udp[1],d=udp[2],pu=udp[3],x=k,r=r,t=t,n=n,amr=T,put=T)

matplot(s0,sol4, type = "l", lwd = 2, lty = 1, main = "European vs American Put",
        xlab="s0", ylab = "put Price", col = c(1:2))
legend("topright", legend = c("European","American"), lwd = 2, col = c(1:2))


#Question 5
trinomEcall=function(u,d,pu,pd,s0,k,r,t,n,log=F){ 
  pm=1-pu-pd
  st=numeric(2*n+1)
  if(log){
    st[1]=exp(log(s0)+u*n)
    for(i in 2:(2*n+1)){st[i]=st[i-1]*exp(d) }  #last period stock prices for log price process
  } else{
    st[1]=s0*u^n
    for(i in 2:(2*n+1)){ st[i]=st[i-1]*d }  #last period stock prices
  }
    c=pmax(st-k,0)   #last period call prices
    for(j in (n-1):0){    #move backwards in steps
      for(i in 1:(2*j+1)){   
        c[i]=(pu*c[i]+pm*c[i+1]+pd*c[i+2])*exp(-r*t/n)  
      }
    }
  return(c[1])
}

q5a <- function(r, sigma, t, n){
  dt <- t/n
  d <- exp(-sigma*sqrt(3*dt))
  u <- 1/d
  pd <- (r*dt*(1-u)+(r*dt)^2+sigma^2*dt)/((u-d)*(1-d))
  pu <- (r*dt*(1-d)+(r*dt)^2+sigma^2*dt)/((u-d)*(u-1))
  return(c(u,d,pu,pd))
}
q5b <- function(r, sigma, t, n){
  dt <- t/n
  d <- -sigma*sqrt(3*dt)
  u <- sigma*sqrt(3*dt)
  pd <- .5*((sigma^2*dt+(r-sigma^2/2)^2*dt^2)/(u^2)-((r-sigma^2/2)*dt)/(u))
  pu <- .5*((sigma^2*dt+(r-sigma^2/2)^2*dt^2)/(u^2)+((r-sigma^2/2)*dt)/(u))
  return(c(u,d,pu,pd))
}

r <- .05
sig <- .24
s0 <- 32
k <- 30
t <- 0.5
n <- c(10,20,40,80,100,200,500)
sol5=matrix(0,length(n),2)
for (i in 1:length(n)){
  udp1 = q5a(r, sig, t, n[i])
  sol5[i,1] = trinomEcall(udp1[1], udp1[2], udp1[3], udp1[4], s0, k, r, t, n[i])
  
  udp2 = q5b(r, sig, t, n[i])
  sol5[i,2] = trinomEcall(udp2[1], udp2[2], udp2[3], udp1[4], s0, k, r, t, n[i],T)
}

matplot(n,sol5, type = "l", lwd = 2, lty = 1, main = "Trinomial Convergence Rates", xlab="steps",
         ylab = "Call Price", col = c(1:2))
abline(h=bscall(s0,k,sig,r,t,0), lty = 5)
legend("topright", legend = paste0("(",letters[1:2],")"), lwd = 2, col = c(1:2))


#Question 6
#from P3
halton = function(n, base){
  seq = rep(0, n)
  numBits = 1+ceiling(log(n)/log(base))
  vetBase = base^(-(1:numBits))
  workVet = rep(0, numBits)
  for(i in 1:n){
    j <- 1
    ok <- 0
    while (ok==0) {
      workVet[j] <- workVet[j]+1
      if (workVet[j]<base){
        ok <-  1
      } else {
        workVet[j] <- 0
        j <- j+1
      }
    }
    seq[i] <- sum(workVet*vetBase)
  }
  return(seq)
}

haltonEcall <- function(s0, k, t, r, sig, n, b1, b2){
  
  h1 <- halton(n/2, b1)
  h2 <- halton(n/2, b2)
  z1 <- sqrt(-2*log(h1))*cos(2*pi*h2)
  z2 <- sqrt(-2*log(h1))*sin(2*pi*h2)
  
  wt <- sqrt(t)*cbind(z1, z2)
  st <- s0*exp((r-sig^2/2)*t + sig*wt)
  c <- exp(-r*t)*mean(pmax(st-k,0))
  return(c)
}



