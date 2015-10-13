# Created by Lacey Westphal, March 6 2015
## This script takes multiple timepoint files ordered by RefPos of boxcox transformed (& centered) IPDS
## and performs an ANOVA on the samples, outputting a .txt file with p_values ordered by RefPos
# file RefPos must be intersected, otherwise this will not work.

##Usage: In terminal: python script.py 

##### DO NOT RUN THIS MULTIPLE TIMES WITHOUT DELETING THE OUTPUT FILE. IT APPENDS!! #####

from scipy import stats
import math
import numpy as np


f1 = '8hr_ordered_int_less1000.txt'
f2 = '16hr_ordered_int_less1000.txt'
f3 = '24hr_ordered_int_less1000.txt'
f4 = '48hr_ordered_int_less1000.txt'
f5 = '72hr_ordered_int_less1000.txt'
f6 = '96hr_ordered_int_less1000.txt'
f7 = '120hr_ordered_int_less1000.txt'
output= "ANOVA_pvals_less1000.txt" 

outfile = open(output, 'a')

# set up an empty list, go through each row of the file and append the IPD if the RefPos
# is the same as the current RefPos. Once it hits a new RefPos, yield the already made list
# to the stats function.

def generate_IPDs(file):
    IPD=[] #set-up empty list of IPD
    current_ref = 1.0 #start current_ref at the beginning
    with open(file, 'r') as f:
        f.next() #skip first str (header)
        for row in f:              
            line = row.split() #split file into rows            
            if current_ref == float(line[1]): #if current_ref is equal to value of RefPos
                IPD.append(float(line[2])) #append to list called IPD                  
            else:    #if not equal, yield the current_ref and IPD list to function 
                yield current_ref, IPD
                del IPD[:] #delete contents of IPD lis
                IPD.append(float(line[2])) #but since we have moved to next line, attach that to new IPD
                current_ref = float(line[1])  #and make sure current_ref value is next RefPos value
    yield current_ref, IPD  #this yields last IPD list when eof is reached        
            
                        
    
def generate_stats(): 
     # create output file
    outfile.write("RefPos"+"\t"+"Pvalue"+"\n")
    tipd1 = generate_IPDs(f1) #create generator object
    tipd2 = generate_IPDs(f2)
    tipd3 = generate_IPDs(f3)
    tipd4 = generate_IPDs(f4)
    tipd5 = generate_IPDs(f5)
    tipd6 = generate_IPDs(f5)
    tipd7 = generate_IPDs(f7) #create generator object (so you can use in for-loop)
    for ref1, ipd1 in tipd1: #this will iterate through generator(control)
        ref2, ipd2 = tipd2.next()
        ref3, ipd3 = tipd3.next()
        ref4, ipd4 = tipd4.next()
        ref5, ipd5 = tipd5.next()
        ref6, ipd6 = tipd6.next()
        ref7, ipd7 = tipd7.next()
         #generate next value from generator (timepoint)        
        if ref1==ref2:
        	if ref2 == ref3:
        		if ref3 == ref4:
        			if ref4 == ref5:
        				if ref5 == ref6:
        					if ref6 == ref7:
        						print ref1
        						ANOVA(ipd1, ipd2, ipd3, ipd4, ipd5, ipd6, ipd7, ref1)
           #if they are equal, pass values to function
        else: 
        	print "your references are not aligned"
            #skip control forward generate_IPDs.next()? 

def ANOVA(ipd1, ipd2, ipd3, ipd4, ipd5, ipd6, ipd7, ref1): #calculates z_stat and pvalue from Wilcoxon test, appends to a file
       f_val, p_val = stats.f_oneway(ipd1, ipd2, ipd3, ipd4, ipd5, ipd6, ipd7)
       outfile.write(str(ref1)+"\t"+str(p_val)+"\n")
       return         

           
generate_stats()         
outfile.close()
