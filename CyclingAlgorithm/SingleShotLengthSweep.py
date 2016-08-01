#%matplotlib inline  
#pylab inline
import random
from random import choice
import copy
import matplotlib.pyplot as plt
import sys
from compiler.ast import flatten
from numpy import *
import numpy

arg1 = float(sys.argv[1])%30 + 1
arg2 = floor(float(sys.argv[1]) / 30.)



#SystemLength = 24
#****

#note, for int(arg2 + 3) this starts with 24, for +2 it starts with 12
NumOfLambdas = 2 ** (int(arg2)+3)
#NumOfLambdas = 13
#****

if (arg2%3==0):
    NumOfLambdas = 16
if (arg2%3==1):
    NumOfLambdas = 32
if (arg2%3==2):
    NumOfLambdas = 64
    








SystemLength = 3 * NumOfLambdas


IsingChain = [0 for i in xrange(SystemLength)]
#IsingChain[0] = 1
#IsingChain[3] = 1

Temperature = .07
#Temperature = .1
Delta = 1.0

CreationRate = abs(1./(1-exp(Delta*1.0/Temperature)))
AnnihilationRate = abs(1./(1-exp(-Delta*1.0/Temperature)))
HoppingRate = .01*Temperature

#CorrectionRate = len(SwapProtocol(SystemLength))





def ReturnExcitationInformation(chain):
    ExLocList = []
    ExPairLocList = []
    EmptyLocList = []
    EmptyPairLocList = []
    RightHoppableLocList = []
    LeftHoppableLocList = []
    for i,c in enumerate(chain):
        if c == 1:
            ExLocList.append(i)
            if chain[(i+1)%len(chain)] == 1:
                ExPairLocList.append(i)
            else:
                RightHoppableLocList.append(i)
        else:
            EmptyLocList.append(i)
            if chain[(i+1)%len(chain)] == 0:
                EmptyPairLocList.append(i)
            else:
                LeftHoppableLocList.append((i)%len(chain))
        
    return ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList



def CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate):
    ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = ReturnExcitationInformation(chain)
    Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate + len(ExPairLocList)*AnnihilationRate + \
    (len(EmptyPairLocList))*CreationRate
    
    
    PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate/Norm
    PAnn = len(ExPairLocList)*AnnihilationRate/Norm
    PCre = (len(chain) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*HoppingRate/Norm
    
    return PHop, PAnn, PCre



def AdvanceTime(chain, StartTime, CreationRate, AnnihilationRate, HoppingRate, sector):
    ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = ReturnExcitationInformation(chain)
    Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate + len(ExPairLocList)*AnnihilationRate + \
    (len(EmptyPairLocList))*CreationRate
    
    PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate/Norm
    PAnn = len(ExPairLocList)*AnnihilationRate/Norm
    PCre = (len(chain) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*HoppingRate/Norm
            
    r = random.random()
    DeltaTau = (-1./Norm)*log(r)
    
    chain, CycleTime, NewRates, NoHops, Proceed, sector = CorrectionProtocol(chain, StartTime, StartTime+DeltaTau, CorrectionRate, \
                                                          PHop, PAnn, PCre, CorrectionSwaps, sector)
    #NewRates = False
    #CycleTime = 0
    #NoHops = True
    
    p = int(floor(len(chain)/2.))
    pl,pr = chain[p],chain[p+1] #previous values of chain
    
    if NewRates == False:
        ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = ReturnExcitationInformation(chain)
        Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate + len(ExPairLocList)*AnnihilationRate + \
        (len(EmptyPairLocList))*CreationRate
        
        PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate/Norm
        PAnn = len(ExPairLocList)*AnnihilationRate/Norm
        PCre = (len(chain) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*HoppingRate/Norm
    
    
        
        #print PHop, PAnn, PCre
        #print r
        #Hopping
        if r < PHop:
            HopSite = choice(RightHoppableLocList + LeftHoppableLocList)
            #chain[HopSite] = 0
            #print RightHoppableLocList
            #print LeftHoppableLocList
            #print chain
            #print ExLocList
            #print "*****"
            chain[HopSite] = 0
            if HopSite in RightHoppableLocList:
                chain[(HopSite+1)%len(chain)] = 1
            else:
                chain[(HopSite+1)%len(chain)] = 0
                chain[(HopSite)] = 1
            #print "Hopping!"
            #print chain
    
        #Annihilating
        if (r >= PHop and r < PHop + PAnn):
            AnnihilateSite = choice(ExPairLocList)
            chain[AnnihilateSite] = 0
            chain[(AnnihilateSite+1)%len(chain)] = 0
            #print "Annihilating!"
            #print chain
    
        #Creating
        if (r >= PHop + PAnn):
            CreateSite = choice(EmptyPairLocList)
            chain[CreateSite] = 1
            chain[(CreateSite+1)%len(chain)] = 1
            #print "Creating!"
            #print chain
        
        
    sector = (sector + CheckSector(IsingChain,p,pl,pr))%2
    
    if NoHops or not(Proceed):
        return chain, DeltaTau, sector
    else:
        return chain, CycleTime, sector

def CheckSector(chain,p,pl,pr):
    increment = 0
    
    if chain[p]!=pl and chain[p+1] != pr:
        increment = 1
    
    
    #print p,pl,pr,"\t",chain[p],chain[p+1],"\t",increment
    #print chain
    
    return increment


#Constructs a list with the indices for conditional swaps in the correction protocol
#Convention is that the value at protocol[i] is CSWAPPED with protocol[(i+1)%length]
def SwapProtocol(length):
    sublength = length/2 - 1
    protocol = []
    for i in xrange(length):
        for j in xrange(sublength):
            for k in xrange(sublength - j):
                protocol.append((i+(j+k))%length)
            for k in xrange(sublength - j):
                protocol.append((i+(sublength-k-1))%length)
                
    return protocol



def SwapProtocol2(length):
    sublength = int(math.ceil(length/2))
    subdomain = int(sublength / 2)
    protocol = []
    
    
    
    for c in xrange(4):
        subprotocol = []
        for i in xrange(subdomain-1):
            subprotocol.append((subdomain*(c+1) + i)%length)
        protocol.append(subprotocol)

        for j in xrange(subdomain-1):
            for k in xrange(j+1):
                protocol.append((subdomain*c + k + (subdomain-1) - (j+1))%length)
            protocol.append(subprotocol)
        
    protocol = flatten(protocol)
    return protocol

def SwapProtocol3(length):
    sublength = int(math.ceil(length/2))
    subdomain = int(sublength / 2)
    protocol = []
    
    
    
    for c in xrange(4):
        subprotocol = []
        for i in xrange(subdomain-1):
            for m in xrange(i+1):
                subprotocol.append((subdomain*(c+1) + i - m)%length)
        protocol.append(subprotocol)

        for j in xrange(subdomain-1):
            for k in xrange(j+1):
                protocol.append((subdomain*c + k + (subdomain-1) - (j+1))%length)
            protocol.append(subprotocol)
        
    protocol = flatten(protocol)
    return protocol

def SwapProtocol4(length):
    sublength = int(math.ceil(length/2))
    subdomain = int(sublength / 4)
    protocol = []
    
    
    
    for c in xrange(8):
        subprotocol = []
        for i in xrange(subdomain-1):
            for m in xrange(i+1):
                subprotocol.append((subdomain*(c+1) + i - m)%length)
        protocol.append(subprotocol)

        for j in xrange(subdomain-1):
            for k in xrange(j+1):
                protocol.append((subdomain*c + k + (subdomain-1) - (j+1))%length)
            protocol.append(subprotocol)
        
    protocol = flatten(protocol)
    return protocol

def SwapProtocol5(length):
    sublength = int(math.ceil(length/2))
    protocol = []

    for j in xrange(sublength):
        if j%2==0:
            protocol.append(2*j)
            protocol.append(2*j+1)
            protocol.append(2*j)
            protocol.append(2*j+1)


    return protocol


def SwapProtocol6(lamb):

    protocol = []
    baselist = [3, 4, 3, 1, 3, 4, 3, 0, 1, 3, 4, 3]
    NumOfLambdas = lamb

    TotalLength = 3*NumOfLambdas

    for k in xrange(NumOfLambdas):
        protocol.append([(c + 3*k)%TotalLength for c in baselist])
    
    protocol = flatten(protocol)
    return protocol






def CSwap(chain,i):
    #print i
    if chain[i]!=chain[(i+1)%len(chain)]:
        inter = chain[i]
        chain[i] = chain[(i+1)%len(chain)]
        chain[(i+1)%len(chain)] = inter
    ####print "Swapping at " + str(i) + "!: ",chain
    return chain
    
def CorrectionProtocol(chain, oldtime, newtime, CorrectionRate, PHop, PAnn, PCre, CorrectionSwaps, sector):
    
    #print "What"
    CycleTime = 0
    #PHop, PAnn, PCre = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
    #print PHop, PAnn, PCre
    ProbabilityHasChanged = False
    #NoChange = True
    
    NumberOfSwaps = len(CorrectionSwaps)
    #Need to calculate where the correction protocol currently is:
    CorrectionPeriod = 1./CorrectionRate
    NumberCompletedCycles, CurrentCycleTime = divmod(oldtime, CorrectionPeriod)
    IndexInCycle = int(floor((CurrentCycleTime / CorrectionPeriod) * NumberOfSwaps))
    Proceed = True
    
    if (oldtime + CorrectionPeriod/NumberOfSwaps) > newtime:
        Proceed = False
        #psuccess = (newtime - oldtime) / (CorrectionPeriod/NumberOfSwaps)
        #if random.random() < psuccess:
        #    Proceed = True
            
    
    while(oldtime+CycleTime < newtime and not(ProbabilityHasChanged) and not(PHop == 0) and Proceed == True):# and not(PAnn > 0)):
        p = int(floor(len(chain)/2.))
        pl,pr = chain[p],chain[p+1] #previous values of chain
        ####print "Timing information: ", CycleTime,"\t",oldtime,"\t", newtime,"\t",(newtime-oldtime)-CycleTime,"\t",CorrectionPeriod/NumberOfSwaps
        #chain = CSwap(chain, CorrectionSwaps[IndexInCycle])        
        #parallel?
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + NumberOfSwaps/4)%NumberOfSwaps])
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 2*NumberOfSwaps/4)%NumberOfSwaps])
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 3*NumberOfSwaps/4)%NumberOfSwaps])
        
        #NumParallel = 0       
        #if (SystemLength/3) % 2 == 0:
        #    NumParallel = SystemLength / 6
        #else:
        NumParallel = int(floor(SystemLength/6.))
        for NumPar in xrange(NumParallel):
            chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + NumPar*24)%NumberOfSwaps])  


 
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 2*NumberOfSwaps/4)%NumberOfSwaps])
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 3*NumberOfSwaps/4)%NumberOfSwaps])
 

        sector = (sector + CheckSector(IsingChain,p,pl,pr))%2 
        
        PHopInter, PAnnInter, PCreInter = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
        #print PHopInter, PAnnInter, PCreInter
        if (PHop != PHopInter or PAnn != PAnnInter or PCre != PCreInter):
            ProbabilityHasChanged = True
        IndexInCycle = (IndexInCycle+1)%NumberOfSwaps
        CycleTime+=CorrectionPeriod/NumberOfSwaps
        
    
    NoHops = (PHop == 0)
    ####print "At end of correction", chain
    #print "Starttime: ",oldtime,"Candidate endtime:",newtime,"Cycle endtime:",oldtime+CycleTime
    #print "New rate equations: ", ProbabilityHasChanged, "Nohops: ", NoHops, "Proceeded?", Proceed
    
    return chain, CycleTime, ProbabilityHasChanged, NoHops, Proceed, sector

def CheckState(chain, sector):
    if sum(chain)==0:
        return 2*sector-1
    else:
        return 0


def ProcessTraj(traj,avgtraj,maxtime):
    #print traj
    #avgtraj[0]+=traj[0][1]
    trajindex = 0
    for i, val in enumerate(avgtraj):
        
        if i < len(traj) and i > 0: #safety first!
            
            while trajindex < len(traj) and traj[trajindex][0] < (1.0*maxtime / len(avgtraj))*i:
                #print "Window: ",(1.0*maxtime / len(avgtraj))*i
                trajindex+=1
            avgtraj[i]+=traj[trajindex-1][1]
    
    return avgtraj


Tests = SwapProtocol2(SystemLength)
#print Tests



#CorrectionRate = AnnihilationRate/len(SwapProtocol(SystemLength))/10.
#CorrectionRate = .001


#****
#CorrectionRate = .15*(1.5)**(-(float(sys.argv[1])-1))
#CorrectionRate = .0012
#CorrectionRate = .0015 * (1.094468)**(-(float(arg1)-1))
CRE = .038085*(1. / SystemLength) - .00002055 #Empirically determined maximum correction rate

CorrectionRate = 4.*CRE * (1.094468)**(-(float(arg1)-1))

#CorrectionRate = 0.00000000000000001

#This should put the maximum rate somewhere around the middle of the 30 jobs (because the 1.09 etc. factor is approx 1/4 when arg1 = 15)


#****

#CorrectionSwaps = SwapProtocol5(SystemLength)
CorrectionSwaps = SwapProtocol6(NumOfLambdas)

#print SwapProtocol3(SystemLength)
#print CorrectionSwaps
zeros = []
ones = []
trajs = [0 for i in xrange(100)]
states = []
fliptimes = []
for i in xrange(1000):
    #print i
    IsingChain = [0 for i in xrange(SystemLength)]
    #IsingChain[0] = 1
    #IsingChain[3] = 1

    counter = 0
    sector = 0
    totalt = 0
    
    traj = []
    state = []
    #while(counter < 1000 and totalt < exp(1./Temperature)):
    FlipHasOccurred = False
    
    while(FlipHasOccurred == False):# and totalt < exp(1./Temperature)):
        #print counter
        curstate = CheckState(IsingChain,sector)
        if curstate == 1:
            FlipHasOccurred = True
            fliptimes.append(totalt)
        #if totalt >= 4*9955416204 and curstate != 1:
        #    FlipHasOccurred = True
        #    fliptimes.append(totalt)
        #    print "Doublelifetime, ending"
        
        
        p = int(floor(len(IsingChain)/2.))
        pl,pr = IsingChain[p],IsingChain[p+1] #previous values of chain
        #print "Before correction", IsingChain
        IsingChain, t, sector = AdvanceTime(IsingChain, totalt, CreationRate, AnnihilationRate, HoppingRate, sector)

        
        totalt+=t

        #print sector
        #print totalt, sector
        '''
        traj.append([totalt, sector])
        if sector == 1:
            ones.append(totalt)
        else:
            zeros.append(totalt)

        state.append([t,curstate])
        '''    
        counter+=1
    with open('/home/daniel.freeman/ToricCodeProject2/longrunswitherror_fixed/run.'+str(SystemLength)+'.'+str(sys.argv[1])+'.dat','a+') as f:
        f.write(str(Temperature)+"\t"+str(NumOfLambdas)+"\t"+str(SystemLength)+"\t"+str(CorrectionRate)+"\t"+str(totalt)+"\t"+str(totalt*totalt)+"\n")
 
    #trajs = ProcessTraj(traj, trajs, 40000.)
    #states.append(state)

#print fliptimes
#print states[0]
#timeweighted = [[s[0]*s[1] for s in ss] for ss in states]
#avgstates = [float(sum(s))/len(s) for s in timeweighted]
    
#print trajs
#print avgstates
#a,b=histogram(ones, linspace(100,20000.,21))
#print "done"

#print Temperature,"\t",HoppingRate,"\t",NumOfLambdas,"\t",SystemLength,"\t",CorrectionRate,"\t",numpy.mean(fliptimes)
print "done"
