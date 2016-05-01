import random
import numpy as np
import sys

#Genealogies are of the form (for example): 
#"(((1:0.01262,5:0.01262):0.11568,4:0.12830):0.14755,(3:0.13660,2:0.13660):0.13925)"
#We want to extract the time til MRCA from these genealogies, as well as left and right clade information -- is this true?


NumSeqs = int(sys.argv[1])
PartitionNum = int(sys.argv[2])
growth = float(sys.argv[3])
Theta = float(sys.argv[4])
NumGeneals = int(sys.argv[5])
Divergence = float(sys.argv[6])

path_to_files = "Divergence_" + "{:.1f}".format(Divergence) + "Theta_" + str(Theta)+ "NumGeneals_" + str(NumGeneals)

#Identify left
def returnPairwiseDiffs(NumSeqs, start,z, LeftNumbers):
	seqs = [z[i].strip(' ').strip('\n')  for i in range(start, start + NumSeqs)]
	seqs = [[digit for digit in ele] for ele in seqs]
	RightNumbers = [ele for ele in range(1, NumSeqs+ 1) if ele not in LeftNumbers]	
	diffs = []
	for i in LeftNumbers:
		for j in RightNumbers:
			diffs.append(sum([max(int(a),int(b))-min(int(a),int(b)) for a, b in zip(seqs[i-1],seqs[j-1])]))
	return diffs

TrueMRCA = []
InferredMRCA = []
InferredMRCA_all = []
fp =open(path_to_files+"/SampleGenealogy.txt", "r")
GenealogyFileLines = fp.readlines()

Gen_num=0
for line_num,line in enumerate(GenealogyFileLines):
	if line[0:4] == "time":
		Tmrca = float(line.split("\t")[1])
		TrueMRCA.append(Tmrca)
	elif line[0:4] == "segs":
		segs = int(line.strip(" ").split(":")[1])
		if segs:  #if there are segregating sites, figure out how to split them into left and right clades
			LeftNumbersFile= open(path_to_files+"/LeftNumbers"+str(Gen_num).zfill(4), "r")
                        LeftNumbers=LeftNumbersFile.readlines()[PartitionNum - 1].split(' ')
			LeftNumbers=[int(ele) for ele in LeftNumbers]
			#returns a list of polymorphisms for each left-clade/right-clade pair
			diffs = returnPairwiseDiffs(NumSeqs, line_num + 2,GenealogyFileLines, LeftNumbers)
			#averages the polymorphisms over each left-clade/right-clade pair
			meandiffs = np.mean(diffs)
			
			InferredMRCA.append(meandiffs)
			InferredMRCA_all.append(diffs)
		else:
			InferredMRCA.append(0.0)
			a = random.choice(range(1, NumSeqs))
			InferredMRCA_all.append([0 for i in range(a * (NumSeqs-a))])
		Gen_num+=1
fp.close()

#We write the true MRCA values, and the values obtained by averaging over inferred left and right clades to file
TrueMRCAfp= open(path_to_files+"/20150202_TrueMRCA", "a")
TrueMRCA = ','.join(map(str, TrueMRCA))
TrueMRCAfp.write(TrueMRCA)
TrueMRCAfp.close()

InferredMRCAfp= open(path_to_files+"/20150202_InferredMRCA_L_and_R_unknown", "a")
InferredMRCA = ','.join(map(str, InferredMRCA))
InferredMRCAfp.write(InferredMRCA)
InferredMRCAfp.close()

#Now write to 20150202_AllInferredSep_L_and_R_unknown resampling genealogies so that there's always the same number regardless of partitioning
InferredMRCA_all_same_size = []
for ele in InferredMRCA_all:
	randomSampleNumbers = (NumSeqs + int(NumSeqs % 2))*(NumSeqs - int(NumSeqs % 2))/4 - len(ele)
	if randomSampleNumbers:
		Complete_with_resamples = ele + list(np.random.choice(ele, size = randomSampleNumbers, replace=True))
		InferredMRCA_all_same_size.append(Complete_with_resamples)
	else:
		InferredMRCA_all_same_size.append(ele)

InferredMRCA_all_same_size_to_write = [' '.join(map(str,ele)) for ele in InferredMRCA_all_same_size]
InferredMRCA_all_same_size_to_write = '\n'.join(map(str, InferredMRCA_all_same_size_to_write))
InferredMRCALearner2 = open(path_to_files+"/20150202_AllInferredSep_L_and_R_unknown", "w")
InferredMRCALearner2.write(InferredMRCA_all_same_size_to_write)
InferredMRCALearner2.close()


