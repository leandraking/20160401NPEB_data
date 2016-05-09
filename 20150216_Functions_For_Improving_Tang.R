

Thetai.B<-function(x,Xi){
  return ((x+1)*length(Xi[Xi==x + 1])/length(Xi[Xi==x]))
}



CloserRobbinsSepThenMerged<-function(x,fintabl){
  return(fintabl[floor(x)+1] + (x-floor(x))*(fintabl[ceiling(x)+ 1] - fintabl[floor(x)+1]))
}

#This is not the right prior for a sample greater than 2!!!  Only use for sample size = =2
postMeanFullModel<-function(xi, theta){
  rate= 1/theta
  p = 1 / (1+rate)
  #p =1/3 #because the rate is equal to 1/.5
  r=1
  return (p*xi + r*p)
}



insideSum<-function(n,i,t){
  numerator<-(2*i-1)*(-1)^i* i *(i-1)/2 *exp(-t*i*(i-1)/2) *prod(n-i+1:n)
  denom<-prod(n:n+i-1)
  return(numerator/denom)
}

priorTMRCAnSeqs = function(n,t){
  helper<-sapply(2:n, function(i) insideSum(n,i,t))
  return(sum(helper))
}


PriorTimesLikelihood<-function(n,mut){
  likelihoodMut<-sapply(seq(0.2,10,.1), function(x) dpois(round(mut),x) * .5*priorTMRCAnSeqs(n,x))
  plot(seq(0.2,10,.1),likelihoodMut)
  return(max(likelihoodMut))
}


