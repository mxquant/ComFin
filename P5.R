#ComFin Project 5 LSMC
rm(list = ls())

firstValueRow= function (x) {   #function to return first non-zero values of each row
  cumSumMat <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[2])
  for (i in 1:(dim(x)[1])) {
    cumSumMat[i, ] <- cumsum(x[i, ])
  }
  cumSumMat2 <- cbind(matrix(0, nrow = dim(x)[1], ncol = 1), cumSumMat[, -(dim(x)[2])])
  outmat <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[2])
  for (i in 1:dim(x)[2]) {
    outmat[, i] <- ifelse(cumSumMat2[, i] > 0, 0, x[, i])
  }
  return(outmat)
}

lsmAmerPut=function(s0=40, x=40, sig=0.2, r=0.06, t=1, n=100000, f="M", k=3, d=0){  #f=Hermite("H"), Laguerre("L"), Monomial("M")

  steps=316*t  #sqrt(N) = 316
  st1 <- matrix(NA, nrow = n/2, ncol = steps)
  st2 <- matrix(NA, nrow = n/2, ncol = steps)
  for (i in 1:(n/2)) {  #generate stock prices
    w=rnorm(steps)
    st1[i, ] = s0 * exp(cumsum((r-d-0.5*sig^2)*(t/steps) + (sig*(sqrt(t/steps))*w)))
    st2[i, ] = s0 * exp(cumsum((r-d-0.5*sig^2)*(t/steps) + (sig*(sqrt(t/steps))*(-w))))
  }
  st=rbind(st1,st2)
  
  X <- ifelse(st < x, st, NA)  #in the money stock price matrix
  CF = matrix(pmax(0, x - st), nrow = n, ncol = steps)  #cash flow (exercise Values) matrix
  Y1 = exp(-r * (t/steps))*CF    #discounted Cash flow, for comparison
  Y2 <- cbind((matrix(NA, nrow = n, ncol = steps - 1)), Y1[, steps])  #another CF matrix, to be updated each period
  CV <- matrix(NA, nrow = n, ncol = steps - 1) #expected contiuation values
  
  if (f == "L"){  #define basis function
    L <- function(x){
      out <- cbind(exp(-x/2),exp(-x/2)*(1-x),exp(-x/2)*(1-2*x+x^2/2),exp(-x/2)*(1-3*x+(3*x^2)/2-x^3/6))
      return(out[,1:k])
    }
  } else if (f == "H"){
    L <- function(x){
      out <- cbind(x^0,2*x,4*x^2-2,8*x^3-12*x)
      return(out[,1:k])
    }
  } else if (f == "M"){
    L <- function(x){
      out =cbind(x^0,x,x^2,x^3)
      return(out[,1:k])
    }
  } else {stop("Incorrect basis function (f)! Please enter 'M','H','L'")}
  
  for (i in (steps - 1):1) {
    if(f=="L"){   #for Laguerre, need to regress L1, then add intercept
      reg = switch(k,
                   NULL,
                   lm(Y2[,i+1] ~ L(X[,i])[,1] + L(X[,i])[,2]),
                   lm(Y2[,i+1] ~ L(X[,i])[,1] + L(X[,i])[,2] + L(X[,i])[,3]),
                   lm(Y2[,i+1] ~ L(X[,i])[,1] + L(X[,i])[,2] + L(X[,i])[,3] + L(X[,i])[,4]))
      CV[,i] = t(matrix(reg$coefficients[-1])) %*% t(L(X[,i])) + reg$coefficients[1]
    } 
    else{   #for Monomials and Hermite, no need to regress L1
      reg = switch(k,  
                   NULL,
                   lm(Y2[,i+1] ~ L(X[,i])[,2]),
                   lm(Y2[,i+1] ~ L(X[,i])[,2] + L(X[,i])[,3]),
                   lm(Y2[,i+1] ~ L(X[,i])[,2] + L(X[,i])[,3] + L(X[,i])[,4]))
      CV[,i] = t(matrix(reg$coefficients)) %*% t(L(X[,i]))
    }
    
    CV[,i] = (ifelse(is.na(CV[, i]), 0, CV[, i]))   #fill NA's with 0
    Y2[,i] = ifelse(CF[,i] > CV[,i], Y1[,i], Y2[,i+1]*exp(-r*(t/steps))) #Update Y2
  }
  
  CV <- ifelse(is.na(CV), 0, CV)
  CVp <- cbind(CV, (matrix(0, nrow = n, ncol = 1)))   #last column CV is 0
  POF <- ifelse(CVp > CF, 0, CF)       #payoff matrix, 1 or 0
  FPOF <- firstValueRow(POF)          #clean up payoff matrix after first 1 in each row
  dFPOF <- matrix(NA, nrow = n, ncol = steps)     #final discounted prices matrix
  for (i in 1:steps) {
    dFPOF[, i] <- FPOF[, i] * exp(-1 * t/steps * r * i)
  }
  price <- mean(rowSums(dFPOF))
  return(price)
}

#a) Laguerre
t=c(0.5,1,2)
k=c(2,3,4)
asol=matrix(0,3,3)
for(i in 1:3){
  for (j in 1:3){
    asol[i,j]=lsmAmerPut(t=t[i],f="L",k=k[j])
    print(i,j)
  }
}
asol

#b) Hermite
bsol=matrix(0,3,3)
for(i in 1:3){
  for (j in 1:3){
    bsol[i,j]=lsmAmerPut(t=t[i],f="H",k=k[j])
    print(c(i,j))
  }
}
bsol

#c)
csol=matrix(0,3,3)
for(i in 1:3){
  for (j in 1:3){
    csol[i,j]=lsmAmerPut(t=t[i],f="M",k=k[j])
    print(c(i,j))
  }
}
csol

