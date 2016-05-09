
args <- commandArgs(trailingOnly = TRUE)

NumSeqs<-as.integer(args[1])
originalMatrix <- read.table(args[2], sep="\t")
NumPartitions <- as.integer(args[3])
growth<-as.numeric(args[4])
NumGenealogies<-as.integer(args[5])
Theta<-as.numeric(args[6])
Divergence<-as.numeric(args[7])

thedir<-paste("./Divergence_",format(Divergence, nsmall=1), "Theta_", format(Theta, nsmall=1), "NumGeneals_", NumGenealogies, "/", sep="")

substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}

write(1,file=paste(thedir,"LeftNumbers",as.character(substrRight(args[2], 4)), sep="" ),append=TRUE)
