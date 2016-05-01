#!/bin/bash

numSeqs=${4}   #2
theta=${2}  #0.5
numGeneals=${3} #2000 
NumPartitions=${5}
NumReps=2 #30 #25
growth=${1}
HalfnumSeqs=$(echo .5*$numSeqs | bc)
HalfnumSeqs=${HalfNumSeqs%.*}
HalfnumSeqs=${6}
divergence=.4
divergence=${7}

flpath=Divergence_${divergence}Theta_${theta}NumGeneals_${numGeneals}

mkdir $flpath

for rep in `seq $NumReps`;
do
        rm $flpath/20150202_AllInferredSep_L_and_R_unknown
        rm $flpath/20150202_InferredMRCA_L_and_R_unknown
        rm $flpath/20150312_Formatted_Genealogies1*
        rm $flpath/20150312_Formatted_Genealogies0*
        rm $flpath/20150312_Formatted_Genealogies*
        rm $flpath/LeftNumbers0*
        rm $flpath/LeftNumbers1*
        rm $flpath/LeftNumbers*
        rm $flpath/20150202_TrueMRCA

	#Create genealogies with admixture using msms.  Save to file SampleGenealogy.txt
#         ../msms/bin/msms -N 10000 -ms $numSeqs $numGeneals  -I 2 ${HalfnumSeqs} ${HalfnumSeqs} 0 -es 0 1 .8 -ej 0 2 3 -ej $(echo .5*${divergence}*${theta} | bc) 1 3  -t $theta -L -T >   $flpath/SampleGenealogy.txt
	../msms/bin/msms -N 10000 -ms $numSeqs $numGeneals  -I 2 ${HalfnumSeqs} ${HalfnumSeqs} 0 -es 0 1 .5 -ej 0 2 3 -ej $(echo .5*${divergence}*${theta} | bc) 1 3  -t $theta -L -T >   $flpath/SampleGenealogy.txt
	
	echo "genealogies created with msms -- saved to SampleGenealogy.txt"
	date	

	#Format genealogies for input into R partitioning program.  Each genealogy in SampleGenealogy.txt gets its own genealogy of the form: path_to_files + "/20150312_Formatted_Genealogies" + str(gen_num)
	python  20150312_FormatSeqs.py $numSeqs ${1} $numGeneals $NumPartitions $theta $divergence
	
	echo "formatted genealogies for input into R partitioning program. Created a file for each genealogy of the form "/20150312_Formatted_Genealogies" + str(gen_num)"
	date

	#For each genealogy, create NumPartitions partitions according to the Tang algorithm.  
	#The idea here is that we don't know what the true left clade right clade partitions are for samples greater than 2, so we might want to average over more than one partitions.  In practice, it doesn't really make so much of a difference. 
	for f in $flpath/20150312_Formatted_Genealogies* ; 
	do
		if [ "${numSeqs}" == "2" ]
			then
				Rscript 20150615TangPartitionNumSeqs2.R $numSeqs $f $NumPartitions $growth $numGeneals $theta $divergence
		fi
		
		if [ "${numSeqs}" != "2" ]
			then 
				Rscript 20150318TangPartition.R $numSeqs $f $NumPartitions $growth $numGeneals $theta $divergence
		fi
	done

	echo "partitions made"
	date

	for partition in `seq $NumPartitions`;
	do
		python 20160331_sep_then_merged_L_R_unknown_Mult_partitions_combined.py $numSeqs $partition $growth $theta $numGeneals $divergence
	done

	echo "files for R analysis made"
	date

	Rscript 20160321IsoReg.R $growth $theta $numGeneals $numSeqs $NumPartitions $divergence		
	
	echo "Completed R analysis using 20160321IsoReg.R.  MSEs calculated."
	echo "done with rep" $rep
	date

done

#rm $flpath/20150202_AllInferredSep_L_and_R_unknown
#rm $flpath/20150202_InferredMRCA_L_and_R_unknown
rm $flpath/20150312_Formatted_Genealogies*
rm $flpath/LeftNumbers*
#rm $flpath/20150202_TrueMRCA
