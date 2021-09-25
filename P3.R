#Computational Finance Project 3

library(plyr) #use raply() instead of a double for loop to add speed
#Q1
#function to simulate X2
simX=function(x0){ 
  n=1000  # <---set number of simulations
  
  sim1x=function(x0){  #function to simulate one X2
    dt=0.001    #<---set dt
    x=rep(x0,2/dt+1)
    for (i in 2:length(x)){
      w=rnorm(1)
      x[i]=x[i-1]+(0.2-0.5*x[i-1])*dt + 2/3*sqrt(dt)*w    #Euler's Scheme
    }
    x[length(x)]  #return X2
  }
  
  raply(n,sim1x(x0),.progress = "text") #apply the function n times to get a distribution of x2
}

#function to simulate Y
simY=function(y0){
  n=1000  #<---set number of simulations
  
  sim1y=function(y0){  #function to simulate one y2 and y3
    dt=0.001    #<---set dt
    t=seq(0,3,dt)
    y=rep(y0,3/dt+1)
    for (i in 2:length(y)){
      z=rnorm(1)
      y[i]=y[i-1]+((2/(1+t[i-1]))*y[i-1] + (1+t[i-1]^3)/3)*dt + (1+t[i-1]^3)/3*sqrt(dt)*z
    }
    c(y[2/dt+1],y[length(y)])  #return Y2 and Y3
  }
  
  raply(n,sim1y(y0),.progress = "text") #apply the function n times to get a distribution of Y
}

x2=simX(1)   #<--input x0
y=simY(0.75) #<--input y0
y2=y[,1]
y3=y[,2]

#output
sum(y2>5)/length(y2)    #p1
mean(sign(x2) * abs(x2)^(1/3))  #e1
mean(y3)   #e2
mean((x2*y2)*(x2>1))  #e3


#Q2
simXY=function(x0){ 
  n=1000  # <---set number of simulations
  
  sim1=function(x0){  
    dt=0.001    #<---set dt
    x=rep(x0,3/dt+1)
    y=rep(1,1/dt+1)
    t=seq(0,3,dt)
    for (i in 2:length(x)){
      w=rnorm(1)
      z=rnorm(1)
      x[i]=x[i-1]+ 0.25*x[i-1]*dt + 1/3*x[i-1]*sqrt(dt)*w - 0.75*x[i-1]*sqrt(dt)*z    #Euler's Scheme
      if (i <= length(y)){
        y[i]=exp(-0.08*t[i-1]+1/3*sqrt(dt)*w+0.75*sqrt(dt)*z)  
      }
    }
    c(x[1/dt+1],x[length(x)],y[length(y)])    #return X1,X3 and Y1
  }
  
  raply(n,sim1(x0),.progress = "text") #apply the function n times to get a distribution of Xt and Yt
}
xy=simXY(1) #<--input x0
x1=xy[,1]
x3=xy[,2]
y1=xy[,3]

#output
mean((1+x3)^(1/3))  #e1
mean(x1*y1)  #e2


#Q3 
#a) Monte Carlo
eurocall=function(s0,t,x,r,sig){
  n=100000   #<---set number of simulations
  w=rnorm(n,0,sqrt(t))
  w2=-w
  s1=s0*exp((r-0.5*sig^2)+sig*w)  #Antithetic method
  s2=s0*exp((r-0.5*sig^2)+sig*w2)
  c=exp(-r*t)*(pmax(s1-x,0) + pmax(s2-x,0))/2
  return(mean(c))
}

#b) Black-Scholes
eurocall_bs <- function(s0,t,x,r,sig){
  d1 <- (log(s0/x)+(r+0.5*sig^2)*t)/(sig*sqrt(t))
  d2 <- d1-sig*sqrt(t)
  bs <- s0*pnorm(d1)-x*exp(-r*t)*pnorm(d2)
  return(bs)
}

#c)
greeks=function(dt,s0,t,x,r,sig){
  D=numeric(length(s0))
  G=numeric(length(s0))
  TH=numeric(length(s0))
  V=numeric(length(s0))
  call=sapply(s0, eurocall_bs,t=t,x=x,r=r,sig=sig)   #initial call prices
  for(i in 1:length(s0)){
    D[i]=(eurocall_bs(s0[i]+dt,t,x,r,sig) - call[i])/dt
    G[i]=(eurocall_bs(s0[i]+dt,t,x,r,sig)+eurocall_bs(s0[i]-dt,t,x,r,sig)-2*call[i])/(dt^2)
    TH[i]=(eurocall_bs(s0[i],t-dt,x,r,sig)-call[i])/dt
    V[i]=(eurocall_bs(s0[i],t,x,r,sig+dt)-call[i])/dt
  }
  plot(s0,D,type = "l", col=i, main = "Delta vs S0")
  plot(s0,G,type = "l", col=i, main = "Gamma vs S0")
  plot(s0,TH,type = "l", col=i, main = "Theta vs S0")
  plot(s0,V,type = "l", col=i, main = "Vega vs S0")
}
greeks(dt=0.004,s0=15:25,t=0.5,x=20,r=0.04,sig=0.25)


#Q4
k=48   #<---set K
t=0.5  #<---set maturity
heston=function(rho,r,s0,v0,sig,a,b){  #simulate 1 path
  dt=t/1000  #<--set dt
  sf = rep(s0, t/dt+1) #full truncation
  vf = rep(v0, t/dt+1)
  sp = rep(s0, t/dt+1) #partial truncation
  vp = rep(v0, t/dt+1)
  sr = rep(s0, t/dt+1) #reflection 
  vr = rep(v0, t/dt+1)
  
  for(i in 2:(t/dt+1)){
    z1=rnorm(1)
    z2=rnorm(1)
    w1=sqrt(dt)*z1
    w2=sqrt(dt)*(rho*z1+sqrt(1-rho^2)*z2)
    
    vf[i]=vf[i-1] + a*(b-max(vf[i-1],0))*dt + sig*sqrt(max(vf[i-1],0))*w2
    sf[i]=sf[i-1] + r*sf[i-1]*dt + sqrt(max(vf[i-1],0))*sf[i-1]*w1
    
    vp[i]=vp[i-1] + a*(b-vp[i-1])*dt + sig*sqrt(max(vp[i-1],0))*w2
    sp[i]=sp[i-1] + r*sp[i-1]*dt + sqrt(max(vp[i-1],0))*sp[i-1]*w1
    
    vr[i]=abs(vr[i-1]) + a*(b-abs(vr[i-1]))*dt + sig*sqrt(abs(vr[i-1]))*w2
    sr[i]=sr[i-1] + r*sr[i-1]*dt + sqrt(abs(vr[i-1]))*sr[i-1]*w1
  }
  c1=max(sf[i]-k,0)*exp(-r*t) #full truncation
  c2=max(sp[i]-k,0)*exp(-r*t) #partial truncation
  c3=max(sr[i]-k,0)*exp(-r*t) #reflection
  return(c(c1,c2,c3))
}

sol=raply(1000,heston(-0.6,0.03,48,0.05,0.42,5.8,0.0625),.progress = "text") #simulate 1000 paths
colnames(sol)=c("C1","C2","C3")
colMeans(sol)  
sol

#Q5
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
#a)
a=cbind(runif(100),runif(100))
#b)
b=cbind(halton(100,2),halton(100,7))
#c)
c=cbind(halton(100,2),halton(100,4))
#d)
plot(a[,1],a[,2],main="Uniform")
plot(b[,1],b[,2],main="Halton Base(2,7)")
plot(c[,1],c[,2],main="Halton Base(2,4)")

#e)
#base(2,7)
xy=cbind(halton(1000,2),halton(1000,7))
x=xy[,1]
y=xy[,2]
ev=matrix(0,1000,1000)
for(i in 1:1000){
  for(j in 1:1000){
    ev[i,j]=exp(-x[i]*y[j])*(sin(6*pi*x[i])+sign(cos(2*pi*y[j]))*abs(cos(2*pi*y[j]))^(1/3))
  }
}
sum(ev)/1000000  #correct way

#base(2,5)
xy=cbind(halton(10000,2),halton(10000,5))
x=xy[,1]
y=xy[,2]
mean(exp(-x*y)*(sin(6*pi*x)+sign(cos(2*pi*y))*abs(cos(2*pi*y))^(1/3)))

#base(2,4)
xy=cbind(halton(10000,2),halton(10000,4))
x=xy[,1]
y=xy[,2]
mean(exp(-x*y)*(sin(6*pi*x)+sign(cos(2*pi*y))*abs(cos(2*pi*y))^(1/3)))
