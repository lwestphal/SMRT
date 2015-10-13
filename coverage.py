#coverage.py

#this script will determine coverage by position, use ordered files
timepoint = "8hr_orderedC_int.txt"

def generate_cov(file):
    IPD=[] #set-up empty list of IPD
    current_ref = 2.0 #start current_ref at the beginning
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

def generate_stats(): # create output file
    tipd = generate_cov(timepoint) #create generator object
     #create generator object (so you can use in for-loop)
    for reft, ipdt in tipd: #this will iterate through generator(control)
        print reft,len(ipdt)   
		
generate_stats()
