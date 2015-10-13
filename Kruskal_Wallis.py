# Created by Lacey Westphal, March 6 2015
## This script takes a control file ordered by RefPos and a timepoint file (also ordered) and performs
## a Kruskal_wallis Test on the two samples, outputting a .txt file with p_values ordered by RefPos

##Usage: enter input and output file names here

##### DO NOT RUN THIS MULTIPLE TIMES WITHOUT DELETING THE OUTPUT FILE. IT APPENDS!! #####

from scipy import stats
import math
from sys import argv 

timepoint1 = "8hr_ordered_intC.txt"
timepoint2  = "16hr_ordered_intC.txt"
timepoint3  = "24hr_ordered_intC.txt"
timepoint4  = "48hr_ordered_intC.txt"
timepoint5  = "72hr_ordered_intC.txt"
timepoint6  = "96hr_ordered_intC.txt"
timepoint7  = "120hr_ordered_intC.txt"
output = "KW_test_allC.txt"
 

outfile = open(output, 'a')

# set up an empty list, go through each row of the file and append the IPD if the RefPos
# is the same as the current RefPos. Once it hits a new RefPos, yield the already made list
# to the stats function.

def generate_IPDs(file):
    IPD=[] #set-up empty list of IPD
    current_ref = 2.0 #start current_ref at the beginning
    with open(file, 'r') as f:
        f.next() #skip first str (header)
        for row in f:              
            line = row.split() #split file into rows            
            if current_ref == float(line[2]): #if current_ref is equal to value of RefPos
                IPD.append(float(line[3])) #append to list called IPD                  
            else:    #if not equal, yield the current_ref and IPD list to function 
                yield current_ref, IPD
                del IPD[:] #delete contents of IPD lis
                IPD.append(float(line[3])) #but since we have moved to next line, attach that to new IPD
                current_ref = float(line[2])  #and make sure current_ref value is next RefPos value
    yield current_ref, IPD  #this yields last IPD list when eof is reached        
            
                        
    
def generate_stats(): 
     # create output file
    outfile.write("RefPos"+"\t"+"Pvalue"+"\n")
    tipd1 = generate_IPDs(timepoint1) #create generator object
    tipd2 = generate_IPDs(timepoint2)
    tipd3 = generate_IPDs(timepoint3)
    tipd4 = generate_IPDs(timepoint4)
    tipd5 = generate_IPDs(timepoint5)
    tipd6 = generate_IPDs(timepoint6)
    tipd7 = generate_IPDs(timepoint7) #create generator object (so you can use in for-loop)
    for ref1, ipd1 in tipd1: #this will iterate through generator(control)
        ref2, ipd2 = tipd2.next()
        ref3, ipd3 = tipd3.next()
        ref4, ipd4 = tipd4.next()
        ref5, ipd5 = tipd5.next()
        ref6, ipd6 = tipd6.next()
        ref7, ipd7 = tipd7.next()
         #generate next value from generator (timepoint)        
        if ref1==ref2: #if they are equal, pass values to Wilcoxon function
            print ref1 
            Kruskal(ipd1,ipd2,ipd3,ipd4,ipd5,ipd6,ipd7, ref1)
        else: 
        	print "your references are not aligned"
            #skip control forward generate_IPDs.next()? 

def Kruskal(ipd1, ipd2, ipd3, ipd4, ipd5, ipd6, ipd7, ref1): #calculates z_stat and pvalue from Wilcoxon test, appends to a file
       stat, p_val = stats.kruskal(ipd1, ipd2, ipd3, ipd4, ipd5, ipd6, ipd7)
       outfile.write(str(ref1)+"\t"+str(p_val)+"\n")
       return         

           
generate_stats()         
outfile.close()
