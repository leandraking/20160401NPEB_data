#!/bin/bash

#SBATCH -n 1
#SBATCH -p serial_requeue
#SBATCH --mem=500
#SBATCH -t 0-23:00 
#SBATCH -J NPEBDivergence

#./20150506_MasterBash.sh $1 $2 $3 $4 $5 $6 $7
#./20150615_MasterBash.sh $1 $2 $3 $4 $5 $6 $7
#./20150801_MasterBash.sh $1 $2 $3 $4 $5 $6 $7

./20160315_MasterBash.sh $1 $2 $3 $4 $5 $6 $7
