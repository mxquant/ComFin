
import numpy as np
from numpy import sqrt, log, sin, cos, exp, mean, repeat, var, mat, floor, flip
from scipy.stats import norm, ncx2 #ncx2 is non central chi-square
from scipy.sparse import diags
import matplotlib.pyplot as pltC
from math import factorial as f
from datetime import date
from pandas import DataFrame
import pandas as pd
import matplotlib.pylab as pylab
import os


#### Q1.a

def zcb(r0, sig, k, r_bar, dt, T, face, n_path):
    
    N =  int(T/dt)
    rates = np.zeros((n_path, N+1))
    rates[:, 0] = r0
    z = np.random.randn(n_path, N+1)
    for i in range(N):
        rates[:, i+1] = rates[:, i] + k*(r_bar - rates[:, i])*dt + sig*sqrt(dt)*z[:, i]
    
    R = dt*np.sum(rates[:, 1:N+1], axis=1)
    discountfactor = mean(np.exp(-R))
    return discountfactor*face

dt = 1/360 
r0 = 0.05
sig = 0.10
k = 0.82
r_bar = 0.05
T = 0.5 
n_path = 1000
face = 1000
q1a = zcb(r0, sig, k, r_bar, dt, T, face, n_path)




###### Q1.b
def couponbond(r0, sig, k, r_bar, dt, T, PMT, face, path):
    vec_ZCB = np.vectorize(zcb)
    couponTimes = np.arange(0.5, T+0.1, 0.5)
    priceOneDollar = vec_ZCB(r0, sig, k, r_bar, dt, couponTimes, 1, path)
    
    payments = np.repeat(PMT, couponTimes.shape[0])
    payments[-1] = PMT+face
    
    price = sum(payments*priceOneDollar)
    
    return price

dt = 1/360
r0 = 0.05
sig = 0.1
k = 0.82
r_bar = 0.05
T = 4 
face = 1000
PMT = 30
path = 1000    
q1b = couponbond(r0, sig, k, r_bar, dt, T, PMT, face, path)



#### Q1.c
def EcallZCB(r0, sig, k, r, dt, t, T, face, strike, path):
    n = int(t/dt)  # steps till call option expiry
    rates = np.zeros((path, n+1))
    rates[:, 0] = r0
    z = np.random.randn(path, n+1)
    for i in range(n):
        rates[:, i+1] = rates[:, i] + k*(r_bar - rates[:, i])*dt + sig*sqrt(dt)*z[:, i]
    
    B=(1-np.exp(-k*(T-t)))/k  #B(t, T)
    A=np.exp( ( r_bar-(sig**2/(2*k**2)) ) * ( B-(T-t) ) - ((sig**2)*(B**2)/(4*k)) )
    ZCBt = face*A*np.exp(-B*rates[:, n])   #P(t, T)
    
    #price of call option
    payOff = ZCBt - strike
    payOff = np.where(payOff>0, payOff, 0)
    
    R = dt*np.sum(rates[:, 1:n], axis=1)
    option_price = mean(np.exp(-R)*payOff)
        
    return option_price

dt = 1/360 
path = 1000
r0 = 0.05
sig = 0.1
k = 0.82
r_bar = 0.05
T = 0.5    #bond maturity 
t = 3/12   #Call option maturity
strike = 980
face = 1000 

q1c = EcallZCB(r0, sig, k, r_bar, dt, t, T, face, strike, path)




##### Q1.d
def coupon(r0, sig, k, r_bar, dt, t, T, PMT, face, path):

    couponTimes = np.arange(0.5, T+0.1, 0.5) #semiannual coupon originally
    couponTimes = couponTimes - t            #option is expiring at time t so coupon payment time will closer by t now
    vec_ZCB = np.vectorize(zcb)
    priceOneDollar = vec_ZCB(r0, sig, k, r_bar, dt, couponTimes, 1, path)
    
    payments = np.repeat(PMT, couponTimes.shape[0])
    payments[-1] = PMT+face
    
    price = sum(payments*priceOneDollar)
    return price


def EcallCoupon(r0, sig, k, r_bar, dt, t, T, face, PMT, strike, path):

    #simulating interest rate till option expiry at t
    n = int(t/dt) #steps till option expiry
    rates = np.zeros((path, n+1))
    rates[:, 0] = r0
    z = np.random.randn(path, n+1)
    for i in range(n):
        rates[:, i+1] = rates[:, i] + k*(r_bar - rates[:, i])*dt + sig*sqrt(dt)*z[:, i]
    
    R = dt*np.sum(rates[:, 1:n], axis=1) #this will be used to discount option payOff
    
    #Now for each path of interest rate till option expiry 
    #Again simulating multiple paths to calculate bond value P(t,T) 
    #price of coupon bond at t expiring at T
    bondprice_t = np.zeros((path, 1))
    for i in range(path):
        bondprice_t[i, 0] =  coupon(rates[i, n], sig, k, r_bar, dt, t, T, PMT, face, 500)
    
    payOff = bondprice_t - strike
    payOff = np.where(payOff > 0, payOff, 0)
    option_price = mean(exp(-R)*payOff) 
    return option_price


n_path = 100    
r0 = 0.05
sig = 0.1
k = 0.82
r_bar = 0.05
dt = 1/360 #time step is one day
T = 4 
face = 1000
PMT = 30
t = 3/12
strike = 980

q1d = EcallCoupon(r0, sig, k, r_bar, dt, t, T, face, PMT, strike, n_path)



### Q1.e
    
r0 = 0.05
sig = 0.1
k = 0.82
r_bar = 0.05
dt = 1/360 #time step is one day
T = 4 
face = 1000
PMT = 30
t = 3/12
K = 980

couponTimes = np.arange(0.5, T+0.1, 0.5)
coupons = np.repeat(PMT, couponTimes.shape[0])
coupons[-1] = face+PMT

# now calculating price of each coupon(equivalent to zero bond) at time t (option expiry time) 

# calculating T-T_i or T-t 
B=(1-np.exp(-k*(couponTimes-t)))/k  #B(t, T)
A=np.exp( ( r_bar-(sig**2/(2*k**2)) ) * ( B-(couponTimes-t) ) - ((sig**2)*(B**2)/(4*k)) )

epsilon = 10**(-10)
rstar = 10**(-10)
pi = A*exp(-B*rstar)
error = sum(coupons*pi)-K

while (error >= epsilon):
    rstar = rstar+0.0001
    pi = A*exp(-B*rstar)
    error = sum(coupons*pi)-K

#new strike price corresponding to each coupon (zero coupon bond for that period)
K_i = A*exp(-B*rstar)

#Now pricing each coupon at time t=0, calculating P(t, Ti)
Bi=(1-np.exp(-k*(couponTimes-0)))/k  #B(t, T)
Ai=np.exp( ( r_bar-(sig**2/(2*k**2)) ) * ( Bi-(couponTimes) ) - ((sig**2)*(B**2)/(4*k)) )
P_t_Ti = Ai*exp(-Bi*r0)


#calculating P(t, T) 
B1 = (1/k)*(1-exp((-k)*(t)))
A1 = exp(((r_bar-((sig**2)/(2*(k**2))))*((B1-(t)))) - ((sig**2)*(B1**2)/(4*k)))
P_t_T = (A1*exp(-B1*r0))

sigp = (sig/k) * (1-exp(-k*(couponTimes-t)))*sqrt( (1-exp(-2*k*t))/(2*k) )

d1 = (1/sigp)*log(P_t_Ti/(K_i*P_t_T))+.5*sigp
d2 = d1 - sigp

c= (P_t_Ti*norm.cdf(d1)) - (K_i * P_t_T * norm.cdf(d2))

q1e = sum(coupons*c)



### 2.a

def ZCB_cir(r0, sig, k, r_bar, dt, T, face, path):
    
    N =  int(T/dt)
    rates = np.zeros((path, N+1))
    rates[:, 0] = r0
    z = np.random.randn(path, N+1)
    for i in range(N):
        rates[:, i+1] = rates[:, i] + k*(r_bar - rates[:, i])*dt + sig*sqrt(rates[:, i])*sqrt(dt)*z[:, i]
    
    R = dt*np.sum(rates[:, 1:N+1], axis=1)
    discountfactor = mean(np.exp(-R))
    return discountfactor*face

def call_ZCB_cir(r0, sig, k, r_bar, dt, t, T, face, strike, path):

    #simulating interest rate till option expiry at t
    n = int(t/dt) #steps till option expiry
    rates = np.zeros((path, n+1))
    rates[:, 0] = r0
    z = np.random.randn(path, n+1)
    for i in range(n):
        rates[:, i+1] = rates[:, i] + k*(r_bar - rates[:, i])*dt + sig*sqrt(rates[:, i])*sqrt(dt)*z[:, i]
    
    R = dt*np.sum(rates[:, 1:n], axis=1) #used to discount option payOff
    
    
    #simulating multiple paths to calculate bond value P(t,T) 
    bondprice_t = np.zeros((path, 1))
    for i in range(path):
        bondprice_t[i, 0] =  ZCB_cir(rates[i, n], sig, k, r_bar, dt, T-t, face, 100)
    
    payOff = bondprice_t - strike
    payOff = np.where(payOff > 0, payOff, 0)
    option_price = mean(exp(-R)*payOff) 
    return option_price

n_path = 1000    
r0 = 0.05
sig = 0.12
k = 0.92
r_bar = 0.055
dt = 1/360 #time step is one day
T = 1 #bond maturity 
face = 1000
t = 0.5
strike = 980

q2a = call_ZCB_cir(r0, sig, k, r_bar, dt, t, T, face, strike, n_path)


#### Q2.b
 
dt = 1/360  
r0 = 0.05
sig = 0.12
k = 0.92
r_bar = 0.055
T = 1 #bond maturity 
face = 1000
t = 0.5 #Option maturity
strike = 980

## interest rate range
rmin = 0
rmax = 0.1
delta_r = 0.001
M = int((rmax-rmin)/delta_r)+1
N = int(t/dt)+1 #time steps between 0 and option maturity t
rates = np.linspace(rmax, rmin, M)

## price of zero bond CIR Explicit way
h1 = sqrt((k**2)+(2*(sig**2)))
h2 = (k+h1)/2
h3 = (2*k*r_bar)/sig**2
A = (  h1* exp(h2*(T-t) )/(h2*(exp(h1*(T-t))-1)+ h1))**h3
B = (exp(h1*(T-t))-1)/(h2*( exp(h1*(T-t)) - 1 ) + h1)
bond_price_at_t = face*(A*exp(-B*rates))
bond_price_at_t = bond_price_at_t.reshape(M, 1)


#construct A matrix
i = np.arange(1, M-1, 1) 
Pu = -0.5*dt * ( sig**2*rates[i]/delta_r**2 + k*(r_bar-rates[i])/delta_r )
Pm = 1 + dt* ( sig**2*rates[i]/delta_r**2 + rates[i])
Pd = -0.5*dt * ( sig**2*rates[i]/delta_r**2 - k*(r_bar-rates[i])/delta_r )
Amat = diags([Pu, Pm, Pd], [0, 1, 2], shape=(M-2, M)).todense()
Amat = np.concatenate((Amat[0], Amat, Amat[M-2-1]), axis=0)
Amat[0, 0:3] = [1, -1, 0]
Amat[M-1, M-3:M] = [0, 1, -1]


F = mat(np.zeros((M, N))) #payoff matrix
F[:, :] = bond_price_at_t-strike  # Call option payoff K-S at maturity
F = np.where(F>0, F, 0)

Bmat = np.zeros((M, 1))        

for i in range(N-2, -1, -1):
    Bmat[:, 0] = F[:, i+1]
    Bmat[0, 0] = rates[0] - rates[1]
    Bmat[M-1, 0] = 0
    mat(F)[:, i] = np.linalg.inv(mat(Amat))*Bmat
    
q2b = np.interp(r0, flip(rates, axis=0), flip(F[:, 0], axis=0))


#### Q2.c

r0 = 0.05
sigma = 0.12
k = 0.92
mean_r = 0.055
dt = 1/360 
face = 1000
t = 0.5 #Option maturity
strike = 980
S = 1 #bond maturity 
t = 0  #current time
T = 0.5 #Option maturity

K = 980/1000 

h1 = sqrt( (k**2)+(2*(sigma**2)) )
h2 = (k+h1)/2
h3 = (2*k*mean_r)/sigma**2

A1 = (h1*exp(h2*(S))/(h2*(exp(h1*(S))-1)+ h1))**h3
B1 = (exp(h1*(S))-1)/(h2*( exp(h1*(S)) - 1 ) + h1)
P_t_S = A1*exp(-B1*r0)

### Price of bond P(t, T); t=0
A2 = (h1*exp(h2*(T))/(h2*(exp(h1*(T))-1)+ h1))**h3
B2 = (exp(h1*(T))-1)/(h2*( exp(h1*(T)) - 1 ) + h1)
P_t_T = A2*exp(-B2*r0)

A_T_S = (h1*exp(h2*(S-T))/(h2*(exp(h1*(S-T))-1)+ h1))**h3
B_T_S = (exp(h1*(S-T))-1)/(h2*( exp(h1*(S-T)) - 1 ) + h1)

### Explicit Valuation formula for call option
theta = sqrt( k**2 + 2*sigma**2 )
phi=(2*theta)/( sigma**2 * (exp(theta*(T-t))-1))
psi=(k+theta)/(sigma**2)

r_star = log(A_T_S/K)/B_T_S


x1 = 2*r_star*(phi+psi+B_T_S)
p1 = (4*k*mean_r/(sigma**2))
q1 = (2*(phi**2)*r0*exp(theta*T))/(phi+psi+B)
chi_1 = ncx2.cdf(x1,p1,q1)

x2 = 2*r_star*(phi+psi)
p2 = p1
q2 = (2*(phi**2)*r0*exp(theta*T))/(phi+psi)
chi_2 = ncx2.cdf(x2,p2,q2)


q2c = face * P_t_S *chi_1 - (strike * chi_2 * P_t_T)



####Q3.a

#function to calculate zero coupon bond price using G2++ method
def ZCB_g2pp(x0, y0, r0, phi, rho, a, b, sig, eta, dt, S, face, path):
    N =  int(S/dt) #till bond expiry
    
    x = np.zeros((path, N+1))
    y = np.zeros((path, N+1))
    rates = np.zeros((path, N+1))
    x[:, 0] = x0
    y[:, 0] = y0
    rates[:, 0] = r0
    
    z = np.random.randn(2*path, N+1)
    z1 = np.split(z, 2)[0]
    z2 = np.split(z, 2)[1]

    
    for i in range(N):
        x[:, i+1] = x[:, i] - a*x[:, i]*dt + sig*sqrt(dt)*z1[:, i]
        y[:, i+1] = y[:, i] - b*y[:, i]*dt + eta*( rho*sqrt(dt)*z1[:, i] + sqrt(1-rho**2)*sqrt(dt)*z2[:, i] )
    
    rates = x + y + phi     
    R = dt*np.sum(rates[:, 1:N+1], axis=1)
    discountfactor = mean(np.exp(-R))
    return face*discountfactor

#function to calculate european put option
def Eput(x0, y0, r0, phi0, rho, a, b, sig, eta, strike, dt, T, S, face, path):
    N =  int(T/dt) # till option expiry
    
    x = np.zeros((path, N+1))
    y = np.zeros((path, N+1))
    rates = np.zeros((path, N+1))
    x[:, 0] = x0
    y[:, 0] = y0
    rates[:, 0] = r0
    
    z = np.random.randn(2*path, N+1)
    z1 = np.split(z, 2)[0]
    z2 = np.split(z, 2)[1]

    
    for i in range(N):
        x[:, i+1] = x[:, i] - a * x[:, i]*dt + sig*sqrt(dt)*z1[:, i]
        y[:, i+1] = y[:, i] - b * y[:, i]*dt + eta * ( rho*sqrt(dt)*z1[:, i] + sqrt(1-rho**2 ) * sqrt(dt) * z2[:, i] )
  
    rates = x + y + phi
    R = dt*np.sum(rates[:, 1:N], axis=1) #this will be used to discount option payOff
    
    #Now for each path of interest rate till option expiry 
    #Again simulating multiple paths to calculate bond value P(T, S) 
    #price of coupon bond at T expiring at S, t is 0
    bondprice_t = np.zeros((path, 1))
    for i in range(path):
        bondprice_t[i, 0] =  ZCB_g2pp(x[i, N], y[i, N], rates[i, N], phi, rho, a, b, sig, eta, dt, S-T, face, 100)
    
    payOff = strike- bondprice_t 
    payOff = np.where(payOff > 0, payOff, 0)
    option_price = mean(exp(-R)*payOff) 
    return option_price


x0 = 0
y0 = 0
phi = 0.03
r0 = 0.03
strike = 950
T = 0.5
S = 1
face = 1000
dt = 1/360
rho = 0.7 #corr
a = 0.1
b = 0.3
sig = 0.03
eta = 0.08

n_path = 1000

q3a = Eput(x0, y0, r0, phi, rho, a, b, sig, eta, strike, dt, T, S, face, n_path)


#####Q3 explicit formula
x0 = 0
y0 = 0
eta = 0.08
strike = 950
phi = 0.03
r0 = 0.03
rho = 0.7
a = 0.1
b = 0.3
sig = 0.03
S = 1
T = 0.5
t = 0
face = 1000

#function to generate V
def ZCB_g2pp_explicit(a,b,sig,eta,T,rho, phi):
    V1 = ((sig**2)/a**2) * (T + (2*exp(-a*T)/a) - (exp(-2*a*T)/(2*a)) - (3/(2*a)))
    V2 = ((eta**2)*(T + (2*exp(-b*T)/b) - (exp(-2*b*T)/(2*b)) - (3/(2*b))))/(b**2)
    V3 = ((2*rho*sig*eta)*(T + ((exp(-a*T) - 1)/a) + ((exp(-b*T) - 1)/b) - ((exp(-(a+b)*T) - 1)/(a+b))))/(a*b)
    V = V1 + V2 + V3
    priceZCB = exp((-phi*T) - (((1- exp(-a*T))/a)*x0) - (((1- exp(-b*T))/b)*y0) + (V/2))
    return (priceZCB)

P_t_T = ZCB_g2pp_explicit(a,b,sig,eta,T,rho, phi)
P_t_S = ZCB_g2pp_explicit(a,b,sig,eta,S,rho, phi)


sig1 = ((sig**2)*((1- exp(-a*(S-T)))**2) * (1 - exp(-2*a*T)))/(2*(a**3))
sig2 = ((eta**2)*((1- exp(-b*(S-T)))**2) * (1 - exp(-2*b*T)))/(2*(b**3)) 
sig3 = ((2*rho*sig*eta)*(1- exp(-a*(S-T)))*(1- exp(-b*(S-T)))*(1- exp(-T*(a+b))))/(a*b*(a+b))
sig_square = sig1 + sig2 + sig3

sig = sqrt(sig_square)
K = strike/1000
d1 = ((log((K*P_t_T)/(P_t_S)))/sig) - (sig/2)
d2 = ((log((K*P_t_T)/(P_t_S)))/sig) + (sig/2)
q3b = -face*P_t_S*norm.cdf(d1) + P_t_T*strike*norm.cdf(d2)

