import random
import numpy as np
import sys
#Genealogy = "(((1:0.01262,5:0.01262):0.11568,4:0.12830):0.14755,(3:0.13660,2:0.13660):0.13925)"

NumSeqs = int(sys.argv[1])
growth = float(sys.argv[2])
NumGeneals = int(sys.argv[3])
NumPartitions = int(sys.argv[4])
Theta = float(sys.argv[5])
Divergence = float(sys.argv[6])

path_to_files = "Divergence_" + str(Divergence) + "Theta_" + str(Theta)+ "NumGeneals_" + str(NumGeneals)

#Identify left
fp =open(path_to_files+ "/SampleGenealogy.txt", "r")
z = fp.readlines()
fp.close()

Gen_num=0
for line_num,line in enumerate(z):
	if line[0:4] == "segs":
		a = int(line.strip(" ").split(":")[1])
		if not a:
			Gen_num+=1
		else:
			formattedGenealogies= open(path_to_files + "/20150312_Formatted_Genealogies" + str(Gen_num).zfill(4), "w")
			for l in range(line_num+2, line_num+2 +NumSeqs):
				zstripped= z[l].strip('\n')
				for ele in zstripped:
					formattedGenealogies.write(str(ele) + '\t')
				formattedGenealogies.write('\n')
				
			Gen_num+=1
			formattedGenealogies.close()
