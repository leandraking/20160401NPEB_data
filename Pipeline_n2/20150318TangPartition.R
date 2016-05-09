args <- commandArgs(trailingOnly = TRUE)

NumSeqs<-as.integer(args[1])
originalMatrix <- read.table(args[2], sep="\t")
NumPartitions <- as.integer(args[3])
growth<-as.numeric(args[4])
NumGeneals<-as.integer(args[5])
Theta<-as.numeric(args[6])
Divergence<-as.numeric(args[7])

#setwd(paste("/n/home08/leking/20150317_NPEB/20150331_SinglePartition/Growth_", format(growth, nsmall=1),sep=""))
#setwd(paste("/n/home08/leking/20150317_NPEB/20150415_better_Automation/Growth_", format(growth, nsmall=1), "NumSeqs_",NumSeqs,"NumGeneals_",NumGeneals,"NumPartitions_",NumPartitions,"Theta_",Theta, sep=""))
#thedir<-paste("/n/home08/leking/20150317_NPEB/20150901_Same_as_Compare_w_Bayes_n8/Growth_", format(growth, nsmall=1), "NumSeqs_",NumSeqs,"NumGeneals_",NumGeneals,"NumPartitions_",NumPartitions,"Theta_",format(Theta,nsmall=1),"Divergence_",format(Divergence, nsmall=1),'/', sep="")

thedir<-paste("./Divergence_",format(Divergence, nsmall=1), "Theta_", format(Theta, nsmall=1), "NumGeneals_", NumGeneals, "/", sep="")

originalMatrix[,ncol(originalMatrix)]<-NULL

#Separating into left and right clades using Tang

InitialSeqs<-1:NumSeqs
#NumSeqs = 5

# seq1<-c(1,0,0,0,0)
# seq2<-c(1,0,0,0,0)
# seq3<-c(1,1,0,0,0)
# seq4<-c(0,0,1,0,0)
# seq5<-c(0,0,1,1,1)

#NumMuts<-length(seq1)

#Tree is ((1,2),3)(4,5)

#originalMatrix<-data.frame(t(matrix(c(seq1,seq2, seq3, seq4, seq5), nrow=NumSeqs, ncol=NumMuts)))

NumMuts= dim(originalMatrix)[2]

#Step 1: compute Hamming matrix
HammingMatrix<-matrix(data=NA, nrow=NumSeqs, ncol=NumSeqs)
for (i in 1:NumSeqs){
  for (j in 1:NumSeqs){
    HammingMatrix[i,j] = sum(abs(originalMatrix[i,]-originalMatrix[j,]))
  }
}

#Step 2: Find a pair i*, j* such that Ki*j* greater than or equal to all Kij.
maximalValues<-which(HammingMatrix==max(HammingMatrix), arr.ind =TRUE)
Kistarjstar<-maximalValues[sample(1:nrow(maximalValues),1),] #when there are more than one possible Kistarjstar
LeftClade<-Kistarjstar[1]
RightClade<-Kistarjstar[2]

#Step 3: Partition sequences into left clade and right clade
#3a: Find sequences that need to be partitioned
RemainingSequencesToBePartitioned<-setdiff(InitialSeqs, c(LeftClade,RightClade))

#Shuffle sequences
RemainingSequencesToBePartitioned<-sample(RemainingSequencesToBePartitioned)

#Place sequences into left and right clades according to smallest mean hamming distance
for (i in RemainingSequencesToBePartitioned){
  hammingLeft<-sapply(LeftClade, function(x) {sum(abs(originalMatrix[x,] - originalMatrix[i,]))})
  MeanLeft<-mean(hammingLeft)
  hammingRight<-sapply(RightClade, function(x) {sum(abs(originalMatrix[x,] - originalMatrix[i,]))})
  MeanRight<-mean(hammingRight)
  if (MeanLeft < MeanRight){
    LeftClade<-c(LeftClade, i)
  }
  else { RightClade<-c(RightClade, i)}
}

#Step 4: Compute the average pairwise nucleotide difference
denomD<-length(LeftClade)*length(RightClade)
numeratorD<-sum(HammingMatrix[LeftClade,RightClade])

D<-numeratorD/denomD

#Step 5: Poisson matrix
M<-dpois(HammingMatrix, D)
M<-upper.tri(M, diag=F) * M

#Step 6: compute the likelihood ratio matrix
maxM<-max(M)
L<-M/maxM

SetToZeroIndices<-which((L<.2 & L>0), arr.ind =TRUE)
L[SetToZeroIndices]<-0

#Step 7: Choose a pair i, j according to the multinomial probability proportional to Lij
choiceOfPair<-sample(matrix(1:(NumSeqs*NumSeqs), ncol=NumSeqs), 1, prob=L)
Kinewjnew<-which(matrix(1:(NumSeqs*NumSeqs), ncol=NumSeqs) == choiceOfPair, arr.ind=T)

#Step 8: Repeat step 3 onwards NumPartition times

for (dddd in 1:NumPartitions){
  LeftClade=Kinewjnew[1]
  RightClade=Kinewjnew[2]
  RemainingSequencesToBePartitioned<-setdiff(InitialSeqs, c(LeftClade,RightClade))
  
  #Shuffle sequences
  RemainingSequencesToBePartitioned<-sample(RemainingSequencesToBePartitioned)
  
  for (i in RemainingSequencesToBePartitioned){
    hammingLeft<-sapply(LeftClade, function(x) {sum(abs(originalMatrix[x,] - originalMatrix[i,]))})
    MeanLeft<-mean(hammingLeft)
    hammingRight<-sapply(RightClade, function(x) {sum(abs(originalMatrix[x,] - originalMatrix[i,]))})
    MeanRight<-mean(hammingRight)
    if (MeanLeft < MeanRight){
      LeftClade<-c(LeftClade, i)
    }
    else { RightClade<-c(RightClade, i)}
  }
  
  #print(sort(LeftClade))
  #print(sort(RightClade))
  
  #write this to file
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  write(LeftClade,file=paste(thedir,"LeftNumbers",as.character(substrRight(args[2], 4)), sep="" ),append=TRUE, ncolumns=length(LeftClade))
}
