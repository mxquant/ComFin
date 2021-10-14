#ComFin P6
rm(list = ls())
gbm=function(s0=100, sig=0.2, r=0.06, t=1, steps=365, n=1000, d=0){  #generate GBM stock price
  st <- matrix(NA, nrow = n, ncol = steps)
  for (i in 1:n) {
    st[i, ] = s0 * exp(cumsum((r-d-0.5*sig^2)*(t/steps) + (sig*(sqrt(t/steps))*rnorm(steps))))
  }
  return(st)
}

#Q1 Fixed Strike Lookback Call and put
s0=98
x=100
r=0.03
t=1
sig=seq(0.12,0.48,0.04)
n=10000

c=numeric(length(sig))  #lookback call
for(i in 1:length(sig)){  
  st=gbm(s0,sig[i],r,t,n=n)  #generate GBM stock price
  c[i]=mean(pmax(apply(st,1,max)-x,0))  #compute option value
}
p=numeric(length(sig))  #lookback put
for(i in 1:length(sig)){  
  st=gbm(s0,sig[i],r,t,n=n)
  p[i]=mean(pmax(-apply(st,1,min)+x,0))
}

plot(sig,c,type="l",col="blue",lwd=2,main="Fixed Strike Lookback Call vs Sigma",ylab="values",xlab="sigma")
plot(sig,p,type="l",col="red",lwd=2,main="Fixed Strike Lookback Put vs Sigma",ylab="values",xlab="sigma")


#Q2
Proj6_2func = function(lambda1=0.2, lambda2=0.4, capT=5){
  #set parameters 
  V0 = 20000
  L0 = 22000
  mu = -0.1
  sigma = 0.2
  gamma = -0.4
  r0 = 0.02
  n = capT*12
  delta = 0.25
  alpha = 0.7
  epsilon = 0.95
  beta = (epsilon-alpha)/capT
  r = (r0+delta*lambda2)/12
  pmt = L0*r/(1-1/(1+r)^n)
  a = pmt/r
  b = pmt/(r*(1+r)^n)
  c = 1+r
  
  nstep = 100
  npath = 1000
  dt = capT/nstep
  V = vector("numeric")
  S = firstJump = 0
  V[1] = V0
  loan_Q = value_Q = value_S = 0
  default_option = default_count = expected_time = 0
  
  for(i in 1:npath){
    payoff = 0
    tau = Q = capT
    # time of the first adverse event S
    S = which(rpois(nstep, lambda2*dt) == 1)[1]
    if(is.na(S)){
      S = nstep+1
    }
    S= S*dt
    
    for(j in 1:nstep){
      t = (j-1)*dt
      Lt = a - b*c^(12*t)
      qt = alpha + beta*t
      # when adverse event happened
      if(S < t){
        value_S = V[j]
        tau = S
        break
      } 
      # when collaertized value less than loan
      if(V[j] <= qt*Lt){
        Q = t
        loan_Q = Lt
        value_Q = V[j]
        tau = Q
        break
      }
      V[j+1] = V[j] + V[j]*(mu*dt + sigma*sqrt(dt)*rnorm(1,0,1) + gamma*rpois(1, lambda1*dt))
    }
    # when tau is less than T, then exercise the option
    if(tau < capT){
      default_count = default_count + 1
      expected_time = expected_time + tau
      if(Q < S){
        payoff = max(0, loan_Q - epsilon*value_Q)
        default_option = default_option + payoff*exp(-r0*Q)
      }
      else{
        payoff = abs(a - b*c^(12*S) - epsilon*value_S)
        default_option = default_option + payoff*exp(-r0*S)
      }
    }
  }
  result = numeric(3)
  result[1] = default_option/npath  #option value
  result[2] = default_count/npath   #default prob
  result[3] = expected_time/default_count   #Et
  return(result)
}
Proj6_2func() #return estimated option value, D prob and Et, with lambda1=0.2, lambda2=0.4, T=5

#graphs
lambda1=seq(0.05,0.4,0.05)
lambda2=seq(0.0,0.8,0.1)
t=seq(3,8,1)

sol1=rep(0,5)  #fixing lambda1 at 0.2
for(i in 1:length(lambda2)){
  for(j in 1:length(t)){
    sol1=rbind(sol1,c(Proj6_2func(0.2,lambda2[i],t[j]),t[j],lambda2[i]))
  }
}
sol1=as.data.frame(sol1[-1,])
colnames(sol1)=c("Option_Value","Default_Prob","Et","T","lambda2")

sol2=rep(0,5)   #fixing lambda2 at 0.4
for(i in 1:length(lambda1)){
  for(j in 1:length(t)){
    sol2=rbind(sol2,c(Proj6_2func(lambda1[i],0.4,t[j]),t[j],lambda1[i]))
  }
}
sol2=as.data.frame(sol2[-1,])
colnames(sol2)=c("Option_Value","Default_Prob","Et","T","lambda1")

#a)
library(ggplot2)
library(gridExtra)
a1=ggplot(sol1,aes(sol1$T,Option_Value,group=factor(lambda2),color=factor(lambda2))) + geom_line() +labs(x="T")
a2=ggplot(sol2,aes(sol2$T,Option_Value,group=factor(lambda1),color=factor(lambda1))) + geom_line() +labs(x="T")
grid.arrange(a1,a2,ncol=2)
#b)
b1=ggplot(sol1,aes(sol1$T,Default_Prob,group=factor(lambda2),color=factor(lambda2))) + geom_line() +labs(x="T")
b2=ggplot(sol2,aes(sol2$T,Default_Prob,group=factor(lambda1),color=factor(lambda1))) + geom_line() +labs(x="T")
grid.arrange(b1,b2,ncol=2)
#c)
c1=ggplot(sol1,aes(sol1$T,Et,group=factor(lambda2),color=factor(lambda2))) + geom_line() +labs(x="T")
c2=ggplot(sol2,aes(sol2$T,Et,group=factor(lambda1),color=factor(lambda1))) + geom_line() +labs(x="T")
grid.arrange(c1,c2,ncol=2)

