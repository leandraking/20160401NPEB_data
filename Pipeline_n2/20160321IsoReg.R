library("Iso")

#Get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
growth <- as.double(args[1]) #this is the exponential growth rate
Theta<-as.numeric(args[2]) #The population mutation rate
MultiplyingTrue<-Theta/.5 #MultiplyingTrue * TrueMRCA = avg(InferredMRCA)
NumGenealogies<-as.numeric(args[3]) #The number of genealogies, or the number of independent loci m
NumSeqs<-as.integer(args[4]) #This is the sample size 
NumPartitions<-as.integer(args[5]) #The number of partitions of the data using Tang's method (we average over these)
Divergence<-as.numeric(args[6]) #The time until divergence

#This is the cutoff for the variance estimate.  If variance above this quantity, we ignore this point
CUTOFF<-.2

#This is the directory where all the simulation files are found.
thedir<-paste("./Divergence_",format(Divergence, nsmall=1), "Theta_", format(Theta, nsmall=1), "NumGeneals_", NumGenealogies, "/", sep="")
source("/n/home08/leking/20150317_NPEB/20150216_Functions_For_Improving_Tang.R")

#The true values of the TMRCA are located here.  These are from the program msms.
TrueMRCA=scan(file=paste(thedir,"20150202_TrueMRCA", sep=""), what=numeric(),sep=",")
TrueMRCA<-TrueMRCA[1:NumGenealogies] ### Why??
TrueMRCA<-TrueMRCA*MultiplyingTrue
#We are going to append these true values to a file so that they are all stored
write(TrueMRCA,file=paste(thedir,"Truen",as.character(NumSeqs),
                   "Growth", as.character(growth), "Theta", Theta ,"NumSeqs",NumGenealogies, sep=""), append=TRUE,
      ncolumns=length(TrueMRCA))

#These are the number of mutations at each site.  It is the average number of mutations between left and right clades for each genealogy
InferredMRCA=scan(file=paste(thedir,"20150202_InferredMRCA_L_and_R_unknown", sep=""), what=numeric(),sep=",")
#InferredMRCA<-InferredMRCA[1:length(InferredMRCA) - 1] ### Why?? Is this averaging over all pairwise differences?
InferredMRCA<-apply(matrix(InferredMRCA, nrow=NumGenealogies), 1, mean) #Average over all pairwise differences
#We are going to append these averages over left-clade/right clade pairs to a file so that they are all stored
write(InferredMRCA ,file=paste(thedir,"Inferredn",as.character(NumSeqs),
                   "Growth", as.character(growth), "Theta", Theta, "NumSeqs",NumGenealogies,sep=""), append=TRUE,
      ncolumns=length(TrueMRCA))

#Rows represent left-clade right-clade polymorphisms for each rep
AllInferredMRCA=read.table(file=paste(thedir,"20150202_AllInferredSep_L_and_R_unknown",sep="")) 

iterations = length(AllInferredMRCA[1,]) #This is equal to `rep` in the 20160315_MasterBash.sh file
maxPt = max(AllInferredMRCA) #The max number of polymorphisms seem between a left-clade right-clade pair for any rep
minPt = min(AllInferredMRCA)

TableForEachIteration = data.frame()

#####################
for (i in 1:iterations) { #for each pair in genealogy (this will be equal to 1 if 2 sequences were used to build genealogy)  ##Is this true ??

  #This is a list of number of polymorphisms for a single left-clade right=clade pair from each genealogy
  Inferreds<-AllInferredMRCA[,i]  
    
  #We get an estimate of the variance for observed number of mutations from minPt to maxPt.  This is stored in 
  mx<-sapply(minPt:maxPt,function(x) length(Inferreds[Inferreds==x]))
  numWeights<-mx[-1]
  denomWeights<-mx[-length(mx)]
  InvVarEstimate<-mapply(function(x,y,z) ((z+1)^2*(x^2/y^2)*(1/x + 1/y))^(-1), numWeights, denomWeights, minPt:(maxPt-1))

  #Write the count and variance results to file
  write(mx, file=paste(thedir, "AllVars_diffWeights", sep=""), append=TRUE, ncolumns=length(mx))
  write(InvVarEstimate, file=paste(thedir, "AllVars_diffWeights", sep=""), append=TRUE, ncolumns=length(InvVarEstimate))
       
  #We only consider points with cutoff greater than CUTOFF
  tabl<-sapply(minPt:maxPt, 
               function(x) ifelse(1/InvVarEstimate[x + 1 - minPt] < CUTOFF & !is.na(InvVarEstimate[x+1-minPt]),
                                  (x+1)*mx[x + 2 - minPt] / mx[x+1-minPt], NA))
  #Print some things
  print(minPt:maxPt)
  print(mx)
  print("counts")
  print(1/InvVarEstimate)
  print("variances")
  print(tabl)
  print("tabl -- all points with variance smaller than cutoff")

  #Run isotonic regression over points with small enough inferred variance
  pav<-pava(tabl[!is.na(tabl)], InvVarEstimate[!is.na(tabl)], decreasing=FALSE, long.out=FALSE, stepfun=FALSE)

  #Add these values back to original table
  tabl[!is.na(tabl)]<-pav
  print(tabl)
  print("tabl after isotonic regression")
                   
                   
  TableForEachIteration<-rbind(TableForEachIteration, tabl)
}

#Here we average all of the tables
fintabl<-colMeans(TableForEachIteration, na.rm=TRUE)

#Make sure that this is equal to 1 (this is only the case if n=2 or partition is completely known)
regmodel=lm(InferredMRCA~TrueMRCA + 0)
write(regmodel$coefficients[1],file=paste(thedir,"lm_n",as.character(NumSeqs),
      "Growth", as.character(growth), "Theta", as.integer(Theta),"NumSeqs", NumGenealogies, sep=""), append=TRUE)

#This would be an estimate of the MSE using Robbins method, by simply rounding the left-clade/right-clade average.
#We will prefer MSE.CloserRobbins
MSE.Robbins<-mean((fintabl[-minPt + round(InferredMRCA)+1]-TrueMRCA)^2, na.rm=TRUE)

#This is the MLE MSE, using just the points that had low enough variance that we used them in our NPEB method
#This ensures that the MLE MSE and the Robbins MSE are comparable
MSE.Helper.Inferred<-sapply(InferredMRCA, function(x) ifelse(is.na(fintabl[-minPt + round(x) + 1]), NA, x))
MSE.MLE<-mean((MSE.Helper.Inferred - TrueMRCA)^2, na.rm=TRUE) #Should be same length as MSE.Robbins vector

CloserRobbinsValues<-sapply(InferredMRCA, function(x) CloserRobbinsSepThenMerged(-minPt + x, fintabl))
MSE.CloserRobbins<-mean((CloserRobbinsValues - TrueMRCA)^2, na.rm=TRUE) # +1 added in function for -minPt+ x
#We're only comparing with the Bayes MLE when n==2 (otherwise we have difficulty with the closed form of the posterior)
if (NumSeqs==2){
  MSE.Bayes<-mean((mapply(postMeanFullModel, MSE.Helper.Inferred, Theta) - TrueMRCA)^2, na.rm=TRUE)
}

#We consider genealogies with InferredMRCA >0 to be able to compare them to GeneTree estimates
MSE.MLE.NonZero<-mean((MSE.Helper.Inferred[InferredMRCA>0] - TrueMRCA[InferredMRCA>0])^2, na.rm=TRUE)
MSE.CloserRobbins.NonZero<-mean((sapply(InferredMRCA[InferredMRCA>0], function(x) CloserRobbinsSepThenMerged(-minPt+x, fintabl)) - TrueMRCA[InferredMRCA>0])^2, na.rm=TRUE)
#We're only comparing with the Bayes MLE when n==2 (otherwise we have difficulty with the closed form of the posterior)
if (NumSeqs==2){
  MSE.Bayes.NonZero<-mean((mapply(postMeanFullModel, MSE.Helper.Inferred[InferredMRCA>0], Theta) - TrueMRCA[InferredMRCA>0])^2, na.rm=TRUE)
}

#We write all the values gathered so far to file L_and_R_unknownSeparateThenMergedData_n... (and print them too!)
ForTheRecord<-c(MSE.CloserRobbins, MSE.MLE, MSE.Bayes, MSE.CloserRobbins.NonZero, MSE.MLE.NonZero, MSE.Bayes.NonZero, sum(!is.na(MSE.Helper.Inferred)), length(unique(InferredMRCA)))
write(ForTheRecord,file=paste(thedir,"L_and_R_unknownSeparateThenMergedData_n",as.character(NumSeqs),
                    "Growth", as.character(growth), "Theta", Theta ,"NumSeqs",NumGenealogies,"_diffWeights", sep=""), append=TRUE,
       ncolumns=length(ForTheRecord))

print(c("MSE.CloserRobbins", "MSE.MLE", "MSE.Bayes", "MSE.CloserRobbins.NonZero", "MSE.MLE.NonZero", "MSE.Bayes.NonZero", "sum(!is.na(MSE.Helper.Inferred))", "length(unique(InferredMRCA))"))
print(ForTheRecord)

#Write all the Closer Robbins values to file CloserRobbinsn...
write(CloserRobbinsValues,file=paste(thedir,"CloserRobbinsn",as.character(NumSeqs),
                   "Growth", as.character(growth), "Theta", as.numeric(Theta),"NumSeqs",NumGenealogies, "_diffWeights", sep=""), append=TRUE,
      ncolumns=length(TrueMRCA))




























