#Feb 17 2015, LW
#This code has been written as a good exercise to teach me how to do stuff in python
#Also, it is supposed to trim down my dataset so that I can actually work with it in R
# This script takes the file from integrate5_2PS.py and only prints out the columns I want
# in this case, that is the ReadNo, RefPos, and IPD. Make sure to > to a new file in terminal
file = open("24hr_integrateT.txt", 'r')
for line in file.readlines():
	line=line.split("\t")
	print line[0]+"\t"+line[2]+"\t"+line[6]
file.close()
