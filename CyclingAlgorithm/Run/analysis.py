import os
import sys
from os import path

import numpy

#directory_path = './templengthsweep1'
#directory_path = './lambdasweep3'
#directory_path = './seriallambdasweep1'
#directory_path = '../longrunswitherror_fixed'
directory_path = './lambdasweep_fixed'
directory_path = sys.argv[1]


files = [x for x in os.listdir(directory_path) if path.isfile(directory_path+os.sep+x)]

print len(files)

for f in files:
        with open(directory_path + '/' + str(f)) as r:

                lines = [line.rstrip('\n') for line in r]
                #if lines[0].split("\t")
                avgtime = [float(l.split("\t")[4]) for l in lines]
                avgtimesqr = [float(l.split("\t")[5]) for l in lines]
                #print numpy.mean(avgtime)
                #print numpy.mean(avgtimesqr)
                error = numpy.sqrt(numpy.mean(avgtimesqr) - numpy.mean(avgtime)**2.0)/numpy.sqrt(len(lines))
                #print lines[0].split("\t")[0] + "\t" + lines[0].split("\t")[1] + "\t" + lines[0].split("\t")[3]+ "\t" + str(numpy.mean(avgtime)/(10**9.))
                print lines[0].split("\t")[0] + "\t" + lines[0].split("\t")[1] + "\t" + lines[0].split("\t")[3]+ "\t" + str(numpy.mean(avgtime)/(2286173128.)) + "\t" + str(error/2286173128.) + "\t" + str(error/numpy.mean(avgtime))

                #print lines[0].split("\t")[0] + "\t" + lines[0].split("\t")[1] + "\t" + lines[0].split("\t")[3]+ "\t" + str(numpy.mean(avgtime)) + "\t" + str(error/numpy.mean(avgtime))

                #print lines[0].split("\t")[0] + "\t" + lines[0].split("\t")[1] + "\t" + lines[0].split("\t")[3]+ "\t" + str(numpy.mean(avgtime)) + "\t" + str(error)

                #'''str((error/numpy.sqrt(len(lines)))/numpy.mean(avgtime)) + "\t" +'''
