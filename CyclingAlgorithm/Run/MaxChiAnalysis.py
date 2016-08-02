import sys
import numpy as np
import os

files = os.listdir('maxlifetimesearch/.')
#print files[0]
os.chdir('maxlifetimesearch')

datadic = {}

lengths = []

for f in files:
    ff = np.loadtxt(f)
    if f.split(".")[1] not in lengths:
        lengths.append(f.split(".")[1])
    #print ff[0]
    chi = ff[0][3]
    lifetimes = [fff[4] for fff in ff]
    squarelifetimes = [fff[5] for fff in ff]
    
    error = (np.sqrt(np.average(squarelifetimes) - (np.average(lifetimes)**2.0))/np.sqrt(len(ff)))
    
    
    #print ff[0:2]
    #print error
    #print error/np.average(lifetimes)
    try:
        datadic[f.split(".")[1]].append([chi,np.average(lifetimes),error])
    except KeyError:
        datadic[f.split(".")[1]]=[]
        datadic[f.split(".")[1]].append([chi,np.average(lifetimes),error])

from scipy.optimize import curve_fit

def func(x, a,b,c):
    return a*np.exp(-b*(x-c)**2.0)

for l in lengths:

    
    #the curve fit function doesn't like fitting with large values
    #in the exponential, so I normalized the data before the fit
    xdat = [d[0] for d in datadic[l]]
    xdat = [x*1000 for x in xdat]
    ydat = [d[1] for d in datadic[l]]
    maxy = max(ydat)
    ydat = [y/maxy for y in ydat]

    popt, pcov = curve_fit(func, xdat, ydat,p0=[1.0,1,1])
    #print popt,pcov
    yfit = [func(x,popt[0],popt[1],popt[2])*maxy for x in xdat]
    ydat = [d[1] for d in datadic[l]]

    maxindex = np.argmax(yfit)
    xdat = [d[0] for d in datadic[l]]
    #the error here is actually an upperbound--it's just the distance between
    #sampled x-points.  In reality, the error is probably a bit lower, because
    #I'm effectively centroiding.
    print xdat[maxindex], yfit[maxindex], xdat[maxindex]-xdat[maxindex-1],l
    #print xdat[np.argmax(ydat)], max(ydat), xdat[maxindex]-xdat[maxindex-1],l
    
#lengths
#for d in datadic["27"]:
#    print str(d[0]) + "\t" + str(d[1])
