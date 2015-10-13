# Created by Lacey Westphal and Peter Sauvey, August 24 2015
## This script randomly samples IPD for a timepoint/RefPos and then feeds those into Wilcoxon
## It will run 100 times and output the "insig" RefPos for each test into a dataframe
##Usage: In terminal: python script.py "control file" "timepoint file" "output file" (must be in order)

##### DO NOT RUN THIS MULTIPLE TIMES WITHOUT DELETING THE OUTPUT FILE. IT APPENDS!! #####

from scipy import stats
import math
from sys import argv 
import numpy
import qvalue
import pandas
from Queue import Queue
from threading import Thread
numpy.set_printoptions(threshold='nan')

timepoint = "/Users/lwestpha/Desktop/Lacey/Bioinformatics/8hr_orderedGATC.txt"
control= "/Users/lwestpha/Desktop/Lacey/Bioinformatics/dam_orderedGATC.txt"


ref = []


def generate_IPDs(file):
    IPD=[] #set-up empty list of IPD
    current_ref = 620.0 
    #start current_ref at the beginning
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
         
def generate_DF(ref):    
    ref_arr = numpy.asarray(ref)
    qval = qvalue.estimate(ref_arr[:,[1]])
    q_val = numpy.asarray(qval)
    df_all  = numpy.append(ref_arr, q_val, 1)
    insig_ref = df_all[df_all[:,2] > 0.2]
    print df_all


def generate_stats(): 
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
        if refc==reft: 
            #if they are equal, pass values to Wilcoxon function
            
            for a in Wilcoxon(ipdc, ipdt):
                ref.append([refc,a]) #this will return a 2d list with a row with the refc and every pval
            
    generate_DF(ref)
    
    
def Wilcoxon(ipdc, ipdt): #calculates z_stat and pvalue from Wilcoxon test, appends to a file
        pv_list = []
        if len(ipdc) > len(ipdt):
        	sample = len(ipdt)
        else:
        	sample = len(ipdc)
        for r in range(100): 
           sampled_timepoint = numpy.random.choice(ipdt, sample, replace=True)
           sampled_control = numpy.random.choice(ipdc, sample, replace=True)
           z_stat, p_val = stats.ranksums(sampled_control, sampled_timepoint)  
           pv_list.append(p_val)
        return pv_list         
           
            
generate_stats()
