import numpy as np
from numpy import sqrt, log, sin, cos, exp, mean, repeat, var, mat, floor, flip
from scipy.stats import norm, ncx2 #ncx2 is non central chi-square
from scipy.sparse import diags
import matplotlib.pyplot as plt
from math import factorial as f
from datetime import date
from pandas import DataFrame
import pandas as pd
import matplotlib.pylab as pylab
from scipy.optimize import minimize
from scipy import optimize


def CIR10(r, k, sig, rbar):
    T = 10
    h1 = np.sqrt((k**2)+(2*(sig**2)))
    h2 = (k+h1)/2
    h3 = (2*k*rbar)/(sig**2)
    A = ((h1*np.exp(h2*T))/((h2*(np.exp(h1*T)-1))+h1))**h3
    B = (np.exp(h1*T)-1)/((h2*(np.exp(h1*T)-1))+h1)
    bond_price = A * np.exp(-B*r)
    r_10 = (-1/T)*(np.log(bond_price))
    return(r_10)

def Numerix_MBS_CIR (WAC, T, notional, r0, sig, kappa, rbar, Npath):
    
    # daily frequency
    freq = 360
    steps = int(T*freq)
    dt = 1/freq

    r_daily = np.zeros((Npath, steps+1))
    r_daily[:, 0] = r0
    Z = np.random.randn(Npath, steps)
    
    for i in range(steps):
        r_daily[:,i+1] = r_daily[:,i] + kappa*(rbar - r_daily[:,i])*dt + sig*np.sqrt(np.abs(r_daily[:,i])) * np.sqrt(dt)*Z[:,i]
    
    month = np.linspace(0, T*360, 361).astype(int)
    r = r_daily[:, month]
    

    # find 10Y UST yield
    r_10 = CIR10(r, kappa, sig, rbar)


    # loop through months 
    freq = 12
    steps = int(T*freq)
    dt = 1/freq
    
    CPR = np.zeros((Npath, steps+1))
    RefinanceIncentive = np.zeros((Npath, steps+1))
    Burnout = np.zeros((Npath, steps+1))
    Seasoning = np.minimum(1, np.linspace(0, steps, steps+1)/30).reshape(1, steps+1)
    Seasonality = [0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98]
    
    c = np.zeros((Npath, steps+1))
    IP = np.zeros((Npath, steps+1))
    PV = np.zeros((Npath, steps+1))
    MP = np.zeros((Npath, steps+1))
    SP = np.zeros((Npath, steps+1))
    PP = np.zeros((Npath, steps+1))
    
    PV[:, 0] = notional
    
    for i in range(steps):
        # CPR
        RefinanceIncentive[:, i+1] = 0.28 + 0.14*np.arctan(-8.57 + 430*(WAC - r_10[:, i]))
        Burnout[:, i+1] = 0.3 + 0.7* (PV[:, i]/PV[:, 0])
        CPR[:, i+1] = RefinanceIncentive[:, i+1] * Burnout[:, i+1] * Seasoning[0, i+1] * Seasonality[np.mod(i+1,12) - 1]
        
        IP[:, i+1] = PV[:, i] * (WAC/12)
        MP[:, i+1] = PV[:, i] * (WAC/12) / (1 - 1/((1+(WAC/12))**(steps-i)))
        SP[:, i+1] = MP[:, i+1] - IP[:, i+1]
        PP[:, i+1] = (PV[:, i] - SP[:, i+1]) * (1 - (1 - CPR[:, i+1])**(1/12))
        PV[:, i+1] =  PV[:, i] - (SP[:, i+1] + PP[:, i+1])
        c[:, i+1] = SP[:, i+1] + PP[:, i+1] +  IP[:, i+1]

    R = np.cumsum(r, axis=1)*dt
    discount = np.exp(-R)
    price = np.mean(np.sum(discount[:, :-1] * c[:, 1:], axis=1))
    
    return price, r, c, CPR

#### Q1a
r0 = 0.078
sig = 0.12
kappa = 0.6
rbar = 0.08
market_price = 110000
WAC = 0.08
T = 30
notional = 100000
Npath = 10000

q1a = Numerix_MBS_CIR (WAC, T, notional, r0, sig, kappa, rbar, Npath)[0]


### Q1b
kappa_range = np.linspace(0.3, 0.9, 7)
rbar_range = np.linspace(0.03, 0.09, 7)
price_range = np.zeros((2, 7))

for i in range(7):
    price_range[0, i] = Numerix_MBS_CIR (WAC, T, notional, r0, sig, kappa_range[i], rbar, Npath)[0]
    price_range[1, i] = Numerix_MBS_CIR (WAC, T, notional, r0, sig, kappa, rbar_range[i], Npath)[0]
    

### Q1b
plt.plot(kappa_range, price_range[0,:], 'r--')
plt.title("MBS price vs k")
plt.show()

### Q1c
plt.plot(rbar_range, price_range[1,:], 'b--')
plt.title("MBS price vs r_bar")
plt.tight_layout()
plt.show()


### Q2

def Price_MBS_Numerix(T, WAC, PV0, r0, k, sig, r_bar, spread=0):
    #inputs= input("Enter parameters") for user prompt 
    path = 10000
    delta_t = 1/12
    mortgage_rate = WAC/12
    N = int(30/delta_t)
    
    #Numerix CPR (Conditional Prepayment Rate)
    SY_array = np.array([0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98]) #seasonality
    SY = np.tile(SY_array, T)
    SY = np.insert(SY, 0, 0) #time t=0 is the begining and end of Jan end is t=1
    BU = np.zeros((path,N+1)) #Burnout 
    rate = np.zeros((path,N+1))
    rate10 = np.zeros((path,N+1))
    SG = np.zeros((1,N+1))
    RI = np.zeros((path,N+1))
    CPR = np.zeros((path,N+1))   
    
    #Cash flow estimation
    PV = np.zeros((path,N+1)) #Present value at different time
    CF = np.zeros((path,N+1)) #Cash flow
    DiscountFactors = np.zeros((path,N+1)) ###name it
    IO = np.zeros((path,N+1))
    PO = np.zeros((path,N+1))
    
    PV[:,0] = PV0
    rate[:,0] = r0
    
    z = np.random.randn(path, N+1)
    
     
    #Using CIR Model simulation for interest rate
    for i in range(rate.shape[1]-1):
        #time t=i+1
        #CPR calculation
        rate[:, i+1] = rate[:, i] + k*(r_bar - abs(rate[:, i]))*delta_t + sig*sqrt(abs(rate[:, i]))*sqrt(delta_t)*z[:, i]
        rate10[:, i+1] = CIR10 ( rate[: ,i], k, sig, r_bar)
        RI[:, i+1] = 0.28 + 0.14*np.arctan(-8.57 + 430*(WAC - rate10[:, i+1]) )
        BU[:, i+1] = 0.3 + 0.7 * ( PV[:, i]/PV[:, 0] )
        SG[0, i+1] = min(1,(i+1)/30)
        CPR[:, i+1] = RI[:, i+1] * BU[:, i+1] * SG[0, i+1] * SY[i+1]
        
        #Cash flow calculation
        CF[:, i+1] = ( PV[:, i]*mortgage_rate) / ( 1-(1+mortgage_rate)**(-N+i) ) + \
                     ( PV[:, i] - PV[:, i] * mortgage_rate * ( 1 / ( 1-(1+mortgage_rate)**(-N+i) )-1) ) * ( 1-(1-CPR[:, i+1])**(1/12) )
        
        IO[: ,i+1] = PV[:, i] * mortgage_rate
        PO[ :,i+1] = CF[:, i+1] - IO[:, i+1]
        PV[:, i+1] = PV[:, i] - PO[: ,i+1]
        DiscountFactors[:, i+1] = exp(-delta_t* ( np.sum(rate[:, 0:(i+2)] + spread, axis=1 )) )
    
    discounted_CashFlow = CF*DiscountFactors #initial price
    discountted_IO = IO*DiscountFactors
    discounted_PO = PO*DiscountFactors
    
    MBS_Price = mean(np.sum(discounted_CashFlow, axis=1))
    IO_Price = mean(np.sum(discountted_IO, axis=1))
    PO_Price = mean(np.sum(discounted_PO, axis=1))
    
    return(MBS_Price, IO_Price, PO_Price)

T = 30 
WAC = 0.08
PV0 = 100000 
r0 = 0.078
k = 0.6
r_bar = 0.08
sig = 0.12
marketPrice = 110000

def OAS_spread(spread):
    MBS_price = Price_MBS_Numerix(T, WAC, PV0, r0, k, sig, r_bar, spread)[0]
    error = abs(MBS_price-marketPrice)
    return error


spread = optimize.newton(OAS_spread, 0, tol=1.48e-06)



### Q3
T = 30 
WAC = 0.08
k = 0.6
r_bar = 0.08
sig = 0.12
spread=-0.0123
PV0 = 100000 
r0 = 0.078
y = 0.0005  

P_plus = Price_MBS_Numerix(T,WAC,PV0,r0,k,sig,r_bar,(spread+y))[0]
P_minus = Price_MBS_Numerix(T,WAC,PV0,r0,k,sig,r_bar,(spread-y))[0]
P0 = Price_MBS_Numerix(T,WAC,PV0,r0,k,sig,r_bar,(spread))[0]
Duration = (P_minus-P_plus)/(2*y*P0)
Convexity = (P_minus -2*P0 + P_plus)/(2*(y**2)*P0)/100


### Q4

T = 30 
WAC = 0.08
PV0 = 100000 
r0 = 0.078
k = 0.6
r_bar = 0.08
sig = 0.12

r_bar = np.arange(0.03, 0.09+0.001, 0.01)
v_Price_MBS_Numerix = np.vectorize(Price_MBS_Numerix)
MBS_Price = v_Price_MBS_Numerix(T,WAC,PV0,r0,k,sig,r_bar)

q4_IO = DataFrame(MBS_Price[1], columns=["IO"], index=r_bar)
q4_PO = DataFrame(MBS_Price[2], columns=["PO"], index=r_bar)


line1 = plt.plot(r_bar, MBS_Price[1], 'r--', label="IO")
line2 = plt.plot(r_bar, MBS_Price[2], 'b', label="PO")
plt.title("Numerix MBS IO/PO tranches Price vs r_bar")
plt.xlabel("r_bar")
plt.ylabel("MBS IO and PO price")
plt.legend()