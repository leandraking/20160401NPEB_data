
findDirectlyAboveValue<-function(x, Xi){
  while(TRUE){
    x = x+1
    if (length(Xi[Xi==x])>0 ){
      return(x)
      break
    }
  }
}

findDirectlyBelowValue<-function(x, Xi){
  while(TRUE){
    x = x-1
    if (length(Xi[Xi==x])>0 ){
      return(x)
      break
    }
    if (x<0){
      break
    }
  }
}


GoodTuring<-function(x,Xi){
  b=1
  if (x == max(Xi)){
    return (2*length(Xi[Xi==x])/(2*x - findDirectlyBelowValue(x, Xi)))
  }
  if (x == min(Xi)){
    return (2*length(Xi[Xi==x])/x)
  }
  return(2*length(Xi[Xi==x])/(findDirectlyAboveValue(x, Xi) - findDirectlyBelowValue(x, Xi)))
}

Smoothedx<-function(x, regmodel){
  return(exp(regmodel$coefficients[1] + regmodel$coefficients[2]*x))
}

Thetai.B.Smoothed<-function(x, regmodel){
  return ((x+1)*Smoothedx(x+1,regmodel)/Smoothedx(x,regmodel))
}


Thetai.B<-function(x,Xi){
  return ((x+1)*length(Xi[Xi==x + 1])/length(Xi[Xi==x]))
}

#When put to 50 performed a little less well.  Needs to be smoothed out sooner
Thetai.B.SmoothedAndNot<-function(x, Xi, regmodel){
  if (length(Xi[Xi==x + 1])<80 | length(Xi[Xi==x])<80){
    return(Thetai.B.Smoothed(x, regmodel))
  }
  else return(Thetai.B(x,Xi))
  
}

CloserRobbins<-function(x, Xi, regmodel){
  return(Thetai.B.SmoothedAndNot(floor(x), Xi, regmodel) + (x - floor(x))*
           (Thetai.B.SmoothedAndNot(ceiling(x), Xi, regmodel) - Thetai.B.SmoothedAndNot(floor(x), Xi, regmodel)))
}

CloserRobbinsSepThenMerged<-function(x,fintabl){
  return(fintabl[floor(x)+1] + (x-floor(x))*(fintabl[ceiling(x)+ 1] - fintabl[floor(x)+1]))
}

#This is actually not the right prior for a sample greater than 2!!!
#What is the prior?
postMeanFullModel<-function(xi, theta){
  rate= 1/theta
  p = 1 / (1+rate)
  #p =1/3 #because the rate is equal to 1/.5
  r=1
  return (p*xi + r*p)
}

############################################

#things for the bayesian prior that you are not using now

#Now I need the bayesian estimates for TMRCA for more than 2 seqs
#http://www.stats.ox.ac.uk/~didelot/popgen/lecture5.pdf
#expoLength<-function(i,t){
#  return(i * (i - 1)/2 * exp(-i*(i-1)*t/2))
#}

#prodj<-function(i,n){
#  if (i == 2){
#    return(prod(sapply(3:n, function(j){j * (j-1) / (j * (j-1) - i * (i-1))})))
#  }
#  if (i==n){
#    return(prod(sapply(2:n-1, function(j){j * (j-1) / (j * (j-1) - i * (i-1))})))
#  }
#  if (i > 2 & i < n) {
#    return(prod(sapply(c(2:(i-1), (i+1):n), function(j){j * (j-1) / (j * (j-1) - i * (i-1))})))
#  }
#}  


insideSum<-function(n,i,t){
  numerator<-(2*i-1)*(-1)^i* i *(i-1)/2 *exp(-t*i*(i-1)/2) *prod(n-i+1:n)
  denom<-prod(n:n+i-1)
  return(numerator/denom)
}

priorTMRCAnSeqs = function(n,t){
  helper<-sapply(2:n, function(i) insideSum(n,i,t))
  return(sum(helper))
}

#dpois(round(unique(InferredMRCA)), lambda = 4)

#plot(seq(0,6,.1), sapply(seq(0,6,.1), function(x) .5*priorTMRCAnSeqs(2,x)))
#par(new=T)
#plot(seq(0,6,.1), sapply(seq(0,6,.1), function(x) dexp(x)))
#max(sapply(seq(0,6,.1), function(x) priorTMRCAnSeqs(10,x)))

PriorTimesLikelihood<-function(n,mut){
  likelihoodMut<-sapply(seq(0.2,10,.1), function(x) dpois(round(mut),x) * .5*priorTMRCAnSeqs(n,x))
  plot(seq(0.2,10,.1),likelihoodMut)
  return(max(likelihoodMut))
}


