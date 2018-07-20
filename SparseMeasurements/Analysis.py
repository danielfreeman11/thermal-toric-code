import os
import numpy
import sys


os.chdir('../Runs/')
files = os.listdir('.')

datadict = {}
for j in [1.0,1.5,2.0]:
	for i in xrange(20):
		datadict[str(7*2*i),j]=[]


chi = float(sys.argv[1])


for ff in files:
	with open(ff) as f: 
		lines = f.readlines()
		for line in lines:
			l = line.split("\t")
			datadict[l[2],float(l[3])].append(float(l[4]))


for i in [2,3,4,5,6,7,8,9,10,11,12,13]:
	print 7*i*2,"\t",numpy.mean(datadict[str(7*2*i),chi]),"\t",numpy.std(datadict[str(7*2*i),chi]),"\t",(numpy.std(datadict[str(7*2*i),chi]/numpy.sqrt(len(datadict[str(7*2*i),chi]))))/numpy.mean(datadict[str(7*2*i),chi]),"\t",len(datadict[str(7*2*i),chi])

#os.chdir('../bin/')
#with open("hist.out","a+") as f:
#	for i in datadict["56",1.0]:
#		f.write(str(i) + "\n")
