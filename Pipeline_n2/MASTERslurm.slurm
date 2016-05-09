#!/bin/bash

for growth in 0.0; 
do
	for theta in 0.5 #1.0; #1.0 2.0 4.0;
	do
		for numGeneals in 8000; # 500 1000 4000; #4000
		do
			for numSeqs in 2; #2 #2; # 2 or 8 ... if you change this change half num seqs too
			do 
				for numPartitions in 1;
				do
					for halfNumSeqs in 1; #1; #couldn't figure out how to do int division in bash... 
					do
						for divergence in 4.0 6.0 8.0; #2.0 4.0;
						do
							sbatch submissionScript.slurm $growth $theta $numGeneals $numSeqs $numPartitions $halfNumSeqs $divergence
					#	./20160315_MasterBash.sh $growth $theta $numGeneals $numSeqs $numPartitions $halfNumSeqs $divergence
						done
					done
				done
			done
		done
	done
done
