# Created by Lacey Westphal and Peter Sauvey, March 6 2015
## This script takes a control file ordered by RefPos and a timepoint file (also ordered) and performs
## a Wilcoxon Rank Sum Test on the two samples, outputting a .txt file with p_values ordered by RefPos

##Usage: In terminal: python script.py "control file" "timepoint file" "output file" (must be in order)
##Most common issue: not changing current_ref 
##### DO NOT RUN THIS MULTIPLE TIMES WITHOUT DELETING THE OUTPUT FILE. IT APPENDS!! #####

from scipy import stats
import math
from sys import argv 
import numpy as np

script, control, timepoint, output =argv
 

outfile = open(output, 'a')

# set up an empty list, go through each row of the file and append the IPD if the RefPos
# is the same as the current RefPos. Once it hits a new RefPos, yield the already made list
# to the stats function.

def generate_IPDs(file):
    IPD=[] #set-up empty list of IPD
    current_ref = 1 #start current_ref at the beginning
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
    tipd = generate_IPDs(timepoint) #create generator object
    cipd = generate_IPDs(control) #create generator object (so you can use in for-loop)
    for refc, ipdc in cipd: #this will iterate through generator(control)
        reft, ipdt = tipd.next() #generate next value from generator (timepoint)
        while reft != refc: #if refpos are not equal
            if refc > reft: #skip to next refpos in timepoint if control ref is higher
                print "skipping %r in %s" %(reft, timepoint)
                reft, ipdt = tipd.next()
                
            else: #skip to next refpos in control if timepoint ref is higher
                print "skipping %r in %s" %(refc, control)
                refc, ipdc = cipd.next()        
        if refc==reft: #if they are equal, pass values to Wilcoxon function
            print refc 
            Wilcoxon(ipdc, ipdt, refc)
    
            #skip control forward generate_IPDs.next()? 

def Wilcoxon(ipdc, ipdt, refc): #calculates z_stat and pvalue from Wilcoxon test, appends to a file
    stat, p_val = stats.ranksums(np.random.choice(ipdc, 15), np.random.choice(ipdt, 15))
    outfile.write(str(refc)+"\t"+str(p_val)+"\n")
    return         

           
generate_stats()         
outfile.close()
