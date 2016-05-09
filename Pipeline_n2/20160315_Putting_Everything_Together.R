thedir<-"/n/home08/leking/20150317_NPEB/20160318_Admixture_diff_proportions"

growths<- c(0.0) #4.0
theta<- c(0.5, 1.0,2.0, 4.0) #c(0.5, 1.0, 2.0, 2.5, 5.0)
numGeneals<-c(250, 500, 1000,  4000)
numSeqs<-c(2) #c(2,8) #c(2, 8)
numPartitions<-c(1) #c( 1, 10)
divergence<-c(2.0,4.0) #c(0.0,1.0,3.0)

p<-1
write(c("g", "t", "m", "n", "d","ASEL_CR", "ASEL_MLE", "ASEL_B", "ASEL_MLE_NONZERO", "ASEL_CR_NONZERO", "ASEL_B_NONZERO", "NUM_GREATER_10", "NUM_UNIQUE", "WEIGHTED_MEAN"),paste(thedir,"/20150810_Putting_Everything_Together_diffWeights2", sep=""), ncolumns=14, append=TRUE)

for (g in growths) {
	for (t in theta) {
		for (m in numGeneals){
			for (n in numSeqs){
				for (d in divergence){
					thefile<-paste(thedir,"/Divergence_",format(d, nsmall=1), "Theta_", format(t, nsmall=1), "NumGeneals_", m, "/L_and_R_unknownSeparateThenMergedData_n",n,"Growth",g,"Theta", as.integer(t), "NumSeqs", m ,"_diffWeights",sep="")
					print(thefile)	
					result<-try(EverythingTable<-read.table(thefile));
					if(class(result)=="try-error") next;

					theVarsfile<-paste(thedir, "/Divergence_",format(d, nsmall=1), "Theta_", format(t, nsmall=1), "NumGeneals_", m,  "/AllVars_diffWeights",sep="")
					VarTable<-readLines(theVarsfile)
					counts<-as.vector(sapply(strsplit(VarTable[1], " "), strtoi))
					counts<-counts[-length(counts)]
					print(counts)
					varValues<-as.vector(sapply(strsplit(VarTable[2], " "), as.numeric))
					print(varValues)

					print(weighted.mean(varValues, counts, na.rm=TRUE)) #need to limit it to the above 10
					
					odds<-seq(from=1, to=length(VarTable), by=2)
					
					print("odds")	
					weightedMean<-c()
					for (i in odds){
						counts<-as.vector(sapply(strsplit(VarTable[i], " "), strtoi))
						counts<-counts[-length(counts)]
						varValues<-as.vector(sapply(strsplit(VarTable[i+1], " "), as.numeric))
						weightedMean<-c(weightedMean,weighted.mean(varValues, counts, na.rm=TRUE)) #need to limit it to the above 10
						
					}
					print(mean(weightedMean))

					z<-c(g,t,m,n,d, mean(EverythingTable$V1), mean(EverythingTable$V2), mean(EverythingTable$V3), mean(EverythingTable$V4), mean(EverythingTable$V5), mean(EverythingTable$V6),mean(EverythingTable$V7), mean(EverythingTable$V8),  mean(weightedMean))
					write(z,paste(thedir,"/20150810_Putting_Everything_Together_diffWeights2", sep=""), ncolumns=length(z), append=TRUE)			
				}
			}
		}
	}
}



#growths<-c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)
#n=8

#ReturnMean<-function(growth){
#  tableForGrowth<-read.table(paste("./Growth_", format(growth, nsmall=1), "/L_and_R_unknownSeparateThenMergedData_n8Growth",growth,"Theta0", sep=""))
#  return(c(growth,mean(tableForGrowth$V1),mean(tableForGrowth$V2), mean(tableForGrowth$V3), mean(tableForGrowth$V4)))
#}


#AllMeans<-sapply(growths, ReturnMean)
#print(AllMeans)

#write(t(AllMeans), file="PuttingEverythingTogether.txt", ncolumns=length(AllMeans[1,]))
#write(t(AllMeans), file="PuttingEverythingTogether.txt")
