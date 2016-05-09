# 20160401NPEB_data

The Pipeline_n2 folder corresponds to our method for a sample of size n=2  We run the program 20160315_MasterBash.sh  using the following arguments: growth theta m n 1 n/2 divergence, which have also been described in the main text

We first create genealogies using the program msms (which must be downloaded beforehand).  The following command for example creates m admixed genealogies:
../msms/bin/msms -N 10000 -ms $numSeqs $numGeneals  -I 2 ${HalfnumSeqs} ${HalfnumSeqs} 0 -es 0 1 .5 -ej 0 2 3 -ej $(echo .5*${divergence}*${theta} | bc) 1 3  -t $theta -L -T >   $flpath/SampleGenealogy.txt

 If we don't know what the true left clade right clade partitions are for samples greater than 2, we might want to average over more than one partitions.  The number of parittions is set by the $numPartitions variable.  In practice, it doesn't really make much of a difference, so we always set numPartitions to 1.
 
 The program 20160315_Putting_Everything_Together.R calculates the average MSE for each set of parameters.
 
-----------------
The Pipeline_n8 folder has the instructions for n=8.

-----------------
 The file 20150421NPEBhumanDataNaturePaper-Copy1.ipynb retrieves the distribution of the number of mutations in the neutral loci dataset from Gronau et al (2011).  The file 20160408NPEB_data_IR-Copy1.ipynb uses our method on this data.
-----------------
 Figure3plot.R produces the plot in figure 3.
 
 
