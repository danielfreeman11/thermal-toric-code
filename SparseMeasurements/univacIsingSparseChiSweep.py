#%matplotlib inline  
#pylab inline
import random
from random import choice
import copy
#import matplotlib.pyplot as plt
import sys
from compiler.ast import flatten
from numpy import *
import numpy

#arg1 = float(sys.argv[1])%30 + 1
#arg2 = floor(float(sys.argv[1]) / 30.) + 1
#arg3 = floor(float(sys.argv[1]) / 60.) + 1


#li = 0

#if float(sys.argv[1]) >=30:
#    li = 1
#if float(sys.argv[1]) >=60:
#    li = 2
#li = int(float(sys.argv[1]) // 30.)


#if float(sys.argv[1]) > 120:
#    HoppingRate = .02*Temperature

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
    
    #chain, CycleTime, NewRates, NoHops, Proceed, sector = CorrectionProtocol(chain, StartTime, StartTime+DeltaTau, CorrectionRate, \
    #                                                      PHop, PAnn, PCre, CorrectionSwaps, sector)
    
    chain, CycleTime, NewRates, NoHops, Proceed, sector = WindowMeasurements(chain, StartTime, StartTime+DeltaTau, CorrectionRate, \
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
    
    
        
        ##print PHop, PAnn, PCre
        ##print r
        #Hopping
        if r < PHop:
            HopSite = choice(RightHoppableLocList + LeftHoppableLocList)
            chain[HopSite] = 0
            if HopSite in RightHoppableLocList:
                chain[(HopSite+1)%len(chain)] = 1
            else:
                chain[(HopSite+1)%len(chain)] = 0
                chain[HopSite] = 1
            ##print "Hopping!"
            ##print chain
    
        #Annihilating
        if (r >= PHop and r < PHop + PAnn):
            AnnihilateSite = choice(ExPairLocList)
            chain[AnnihilateSite] = 0
            chain[(AnnihilateSite+1)%len(chain)] = 0
            ##print "Annihilating!"
            ##print chain
    
        #Creating
        if (r >= PHop + PAnn):
            CreateSite = choice(EmptyPairLocList)
            chain[CreateSite] = 1
            chain[(CreateSite+1)%len(chain)] = 1
            ##print "Creating!"
            ##print chain
        
        
    sector = (sector + CheckSector(IsingChain,p,pl,pr))%2
    
    if NoHops or not(Proceed):
        return chain, DeltaTau, sector
    else:
        return chain, CycleTime, sector

def CheckSector(chain,p,pl,pr):
    increment = 0
    
    if chain[p]!=pl and chain[p+1] != pr:
        increment = 1
    
    
    ##print p,pl,pr,"\t",chain[p],chain[p+1],"\t",increment
    ##print chain
    
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

def SwapProtocol7(length, lamb):
    protocol = []
    if length%lamb!=0 or (length/lamb)%2!=0:
        #print "Bad! length lambda mismatch! Must divide the length and quotient must be even!"
        lamb = lamb
    else:
        numofdomains = length / lamb
        for d in xrange(numofdomains):
            for k in xrange(lamb):
                for m in xrange(k):
                    protocol.append((lamb-1-k+m+d*lamb)%length)
                for i in xrange(lamb):
                        for j in xrange(i):
                            protocol.append((lamb+i-j-1+d*lamb)%length)
                

    return protocol    






def CSwap(chain,i):
    ##print i
    if chain[i]!=chain[(i+1)%len(chain)]:
        inter = chain[i]
        chain[i] = chain[(i+1)%len(chain)]
        chain[(i+1)%len(chain)] = inter
    #####print "Swapping at " + str(i) + "!: ",chain
    return chain

def WindowMeasurements(chain, oldtime, newtime, CorrectionRate, PHop, PAnn, PCre, CorrectionSwaps, sector):
    
    #For codes of size (3 + 2n)*k for n>=0,k>=2
    #For now, I'm going to assume n=2
    #So SystemSize = 14,28,56, etc.
    
    
    CycleTime = 0
    ProbabilityHasChanged = False

    
    #NumberOfSwaps = len(CorrectionSwaps)
    NumberOfSwaps = 1
    
    #Need to calculate where the correction protocol currently is:
    CorrectionPeriod = 1./CorrectionRate
    
    #NumberCompletedCycles, CurrentCycleTime = divmod(oldtime, CorrectionPeriod)
    #IndexInCycle = int(floor((CurrentCycleTime / CorrectionPeriod) * NumberOfSwaps))
    IndexInCycle = 0
    
    Proceed = True
    
    if (oldtime + CorrectionPeriod/NumberOfSwaps) > newtime:
        Proceed = False
        #psuccess = (newtime - oldtime) / (CorrectionPeriod/NumberOfSwaps)
        #if random.random() < psuccess:
        #    Proceed = True
    counter2 = 0
    '''#print "Timing check: " + str(oldtime+CycleTime)
    #print "NewTime: " + str(newtime)
    #print "oldtime + corrper" + str(oldtime + CorrectionPeriod/NumberOfSwaps)
    #print oldtime+CycleTime < newtime
    #print not(ProbabilityHasChanged)
    #print not(PHop == 0)
    #print Proceed == True'''

    while(oldtime+CycleTime < newtime and not(ProbabilityHasChanged) and not(PHop == 0) and Proceed == True):# \
        #and counter2 < 100:# and not(PAnn > 0)):
        counter2 += 1
        #print "Cur time: " + str(oldtime + CycleTime) + ", counter:" + str(counter2)
        #print "Newtime: " + str(newtime)
        #print "Chainstate: " + str(chain)
        #print "Sector: " + str(sector)
        p = int(floor(len(chain)/2.))
        pl,pr = chain[p],chain[p+1] #previous values of chain
        
        
        #Check windows, perform simple corrections and simple swaps.
        
        for k in xrange(NumOfLambdas):
            #center defects on measurement sites
            if chain[(7*k-3)%SystemLength]+chain[(7*k-2)%SystemLength]+chain[(7*k-1)%SystemLength]==1:
                chain[(7*k-3)%SystemLength]=0
                chain[(7*k-2)%SystemLength]=1
                chain[(7*k-1)%SystemLength]=0
                
                #print "Centering!"
                
                
            
            #fix simple errors
            if chain[(7*k-3)%SystemLength]+chain[(7*k-2)%SystemLength]+chain[(7*k-1)%SystemLength]==2:
                chain[(7*k-3)%SystemLength]=0
                chain[(7*k-2)%SystemLength]=0
                chain[(7*k-1)%SystemLength]=0
                
                #print "Fixing simple error!"
        
            #for k in xrange(NumOfLambdas):
            if chain[(7*k-2)%SystemLength]==1 and DefectAge[k]==0:
                DefectAge[k]=oldtime+CycleTime
            if chain[(7*k-2)%SystemLength]==0:
                DefectAge[k]=0

                
                
        #Heuristic for error correction: Once Age/Distance exceeds a threshold, pair defects
        #pair defects that are sufficiently old/close:
        
        MaxPairDefects = []
        for ii,d1 in enumerate(DefectAge):
            for jj,d2 in enumerate(DefectAge):
                if d1<d2 and d1*d2!=0:
                    d1i = (7*ii-2)%SystemLength
                    d2i = (7*jj-2)%SystemLength
                    #This calculates the correct pairwise distance accounting for periodic B.C.s
                    errdist = min(abs(d1i-d2i),\
                                  abs(((SystemLength/2.)+d1i)%SystemLength - \
                                      ((SystemLength/2.)+d2i)%SystemLength))
                    MaxPairDefects.append((d1,ii,d2,jj,errdist,d1i,d2i))
        
        MaxPairDefects=sorted(MaxPairDefects, key=lambda x: (-1.0*((errdist*4./7.)**2.0)/\
                              (40.*HoppingRate*abs(oldtime + CycleTime - min(x[0],x[2])))), reverse=True)
        #MaxPairDefects = [(i[0],i[1]) for i in MaxPairDefects]
        for m in MaxPairDefects:
            d1 = m[0]
            i1 = m[1]
            d2 = m[2]
            i2 = m[3]
            
            errdist = m[4]
            
            d1i = m[5]
            d2i = m[6]
            
            #if (d1+d2)*1.0*numpy.exp(-1.*errdist/(3.)) > 4.:
            #print "Fixing?? " + str(.5*abs(d1-d2)-errdist)
            ##print "Agediff: " + str(abs(d1-d2))
            ##print "errdist: " + str(errdist)
            ##print "errdist/hoprate: " + str(errdist/HoppingRate)
            if numpy.random.rand() < min(1, numpy.exp((-1.0*((errdist*4./7.)**2.0)/\
                              (40.*HoppingRate*abs(oldtime + CycleTime - min(d1,d2)))))):
                c1 = 0
                c2 = 0
                #print "Fixing errors!"
                #print "Exponent:  " + str(3.*abs(d1-d2)-errdist)
                #print "Errordist: " + str(errdist)
                #print "Error state before fixing:" + str(chain)
                #print "Fixing index" + str(d1i) + "," + str(d2i)
                #print "Agestate: " + str(DefectAge)
                chain[d1i]=0
                chain[d2i]=0

                DefectAge[((d1i+2)%SystemLength)/7]=0
                DefectAge[((d2i+2)%SystemLength)/7]=0
                
                c1 = d1i
                c2 = d2i
                #print "Error locations: " + str(c1) + "," + str(c2)
                #print "Error min: " + str(min(c1,c2))
                #print "Error max: " + str(max(c1,c2))
                #print "Error disp: " + str(abs(c1-c2))
                #print "Systemdiv: " + str(SystemLength/2)
                if (min(c1,c2) <= SystemLength/2 and max(c1,c2) > SystemLength/2) \
                and abs(c1-c2) <= SystemLength/2:
                    #print "updating sector!"
                    sector = (sector + 1)%2
                #print "Error state after fixing:" + str(chain)
        
                        
                        
                            
        
        
        
        
        
        #NumParallel = NumOfLambdas / 2
        #for NumPar in xrange(NumParallel):
        #    AboutToSwap = CorrectionSwaps[(IndexInCycle + NumPar*2*BaseCycleLength)%NumberOfSwaps]
        #    
        #    #This if statement checks to see if the protocol is about to swap a defect off one of the
        #    #measurement rails.  If it is, it doesn't do the swap (i.e., defects should accumulate on
        #    #the measurement rails)
        #    if (not ((AboutToSwap in Attractors and chain[AboutToSwap] == 1) or \
        #          (AboutToSwap+1)%SystemLength in Attractors and chain[(AboutToSwap+1)%SystemLength] == 1)):
        #        chain = CSwap(chain, AboutToSwap)  
                

                

        #sector = (sector + CheckSector(IsingChain,p,pl,pr))%2 
        
        PHopInter, PAnnInter, PCreInter = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
        ##print PHopInter, PAnnInter, PCreInter
        if (PHop != PHopInter or PAnn != PAnnInter or PCre != PCreInter):
            ProbabilityHasChanged = True
        IndexInCycle = (IndexInCycle+1)%NumberOfSwaps
        CycleTime+=CorrectionPeriod/NumberOfSwaps
        
    
    NoHops = (PHop == 0)
    #####print "At end of correction", chain
    ##print "Starttime: ",oldtime,"Candidate endtime:",newtime,"Cycle endtime:",oldtime+CycleTime
    ##print "New rate equations: ", ProbabilityHasChanged, "Nohops: ", NoHops, "Proceeded?", Proceed
    
    return chain, CycleTime, ProbabilityHasChanged, NoHops, Proceed, sector

    
def SparseMeasurements(chain, oldtime, newtime, CorrectionRate, PHop, PAnn, PCre, CorrectionSwaps, sector):
    ##print "What"
    CycleTime = 0
    #PHop, PAnn, PCre = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
    ##print PHop, PAnn, PCre
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
    counter2 = 0
    
    while(oldtime+CycleTime < newtime and not(ProbabilityHasChanged) and not(PHop == 0) and Proceed == True):# \
        #and counter2 < 100:# and not(PAnn > 0)):
        counter2 += 1
        #print "Cur time: " + str(oldtime + CycleTime) + ", counter:" + str(counter2)
        #print "Newtime: " + str(newtime)
        #print "Chainstate: " + str(chain)
        #print "Sector: " + str(sector)
        p = int(floor(len(chain)/2.))
        pl,pr = chain[p],chain[p+1] #previous values of chain
        
        
        #age Measurements:
        
        for k in xrange(len(Attractors)/2):
            if (chain[Attractors[2*k]]==1 or chain[Attractors[2*k+1]]==1) and DefectAge[k]==0:
                #DefectAge[k]=min(DefectAge[k]+1,SystemLength/2)
                DefectAge[k]=oldtime+CycleTime
            if (chain[Attractors[2*k]]==0 and chain[Attractors[2*k+1]]==0):
                DefectAge[k]=0
            #else:
            #    DefectAge[k]=0
                
                
        #Heuristic for error correction: Once Age/Distance exceeds a threshold, pair defects
        #pair defects that are sufficiently old/close:
        
        MaxPairDefects = []
        for ii,d1 in enumerate(DefectAge):
            for jj,d2 in enumerate(DefectAge):
                if d1<d2 and d1*d2!=0:
                    d1i = Attractors[2*ii]
                    d2i = Attractors[2*jj]
                    #This calculates the correct pairwise distance accounting for periodic B.C.s
                    errdist = min(abs(d1i-d2i),\
                                  abs(((SystemLength/2.)+d1i)%SystemLength - \
                                      ((SystemLength/2.)+d2i)%SystemLength))
                    MaxPairDefects.append((d1,ii,d2,jj,errdist))
        MaxPairDefects=sorted(MaxPairDefects, key=lambda x: (.05*abs(x[0]-x[2])-errdist), reverse=True)
        #MaxPairDefects = [(i[0],i[1]) for i in MaxPairDefects]
        for m in MaxPairDefects:
            d1 = m[0]
            i1 = m[1]
            d2 = m[2]
            i2 = m[3]
            errdist = m[4]
            
            #if (d1+d2)*1.0*numpy.exp(-1.*errdist/(3.)) > 4.:
            ##print "Fixing?? " + str(abs(d1-d2)-errdist/(HoppingRate))
            ##print "Agediff: " + str(abs(d1-d2))
            ##print "errdist: " + str(errdist)
            ##print "errdist/hoprate: " + str(errdist/HoppingRate)
            if numpy.random.rand() < min(1, numpy.exp(.05*abs(d1-d2)-errdist)):
                c1 = 0
                c2 = 0
                #print "Fixing errors!"
                #print "Exponent:  " + str(.05*abs(d1-d2)-errdist)
                #print "Errordist: " + str(errdist)
                #print "Error state before fixing:" + str(chain)
                #print "Fixing index" + str(Attractors[2*i1]) + "," + str(Attractors[2*i2])
                #print "Attractors: " + str(Attractors)
                #print "Agestate: " + str(DefectAge)
                if chain[Attractors[2*i1]]==1:
                    chain[Attractors[2*i1]]=0
                    c1 = Attractors[2*i1]
                else:
                    if chain[Attractors[2*i1 + 1]]==1:
                        chain[Attractors[2*i1 + 1]]=0
                        c1 = Attractors[2*i1 + 1]
                if chain[Attractors[2*i2]]==1:
                    chain[Attractors[2*i2]]=0
                    c2 = Attractors[2*i2]
                else:
                    if chain[Attractors[2*i2 + 1]]==1:
                        chain[Attractors[2*i2 + 1]]=0
                        c2 = Attractors[2*i2 + 1]
                #print "Error locations: " + str(c1) + "," + str(c2)
                #print "Error min: " + str(min(c1,c2))
                #print "Error max: " + str(max(c1,c2))
                #print "Error disp: " + str(abs(c1-c2))
                #print "Systemdiv: " + str(SystemLength/2)
                if (min(c1,c2) <= SystemLength/2 and max(c1,c2) > SystemLength/2) \
                and abs(c1-c2) <= SystemLength/2:
                    #print "updating sector!"
                    sector = (sector + 1)%2
                #print "Error state after fixing:" + str(chain)
                        
            
            
        
        '''for i1,d1 in enumerate(DefectAge):
            for i2,d2 in enumerate(DefectAge):
                if d1<d2 and d1*d2!=0:
                    d1i = Attractors[2*i1]
                    d2i = Attractors[2*i2]
                    #This calculates the correct pairwise distance accounting for periodic B.C.s
                    errdist = min(abs(d1i-d2i),\
                                  abs(d1i - ((SystemLength/2.)+d2i)%SystemLength))
                    
                    
                    if (d1+d2)*1.0/errdist > 10.:
                        c1 = 0
                        c2 = 0
                        #print "Fixing errors!"
                        #print "Error state before fixing:" + str(chain)
                        #print "Fixing index" + str(Attractors[2*i1]) + "," + str(Attractors[2*i2])
                        #print "Attractors: " + str(Attractors)
                        #print "Agestate: " + str(DefectAge)
                        if chain[Attractors[2*i1]]==1:
                            chain[Attractors[2*i1]]=0
                            c1 = Attractors[2*i1]
                        else:
                            if chain[Attractors[2*i1 + 1]]==1:
                                chain[Attractors[2*i1 + 1]]=0
                                c1 = Attractors[2*i1 + 1]
                        if chain[Attractors[2*i2]]==1:
                            chain[Attractors[2*i2]]=0
                            c2 = Attractors[2*i2]
                        else:
                            if chain[Attractors[2*i2 + 1]]==1:
                                chain[Attractors[2*i2 + 1]]=0
                                c2 = Attractors[2*i2 + 1]
                        #print "Error locations: " + str(c1) + "," + str(c2)
                        #print "Error min: " + str(min(c1,c2))
                        #print "Error max: " + str(max(c1,c2))
                        #print "Error disp: " + str(abs(c1-c2))
                        #print "Systemdiv: " + str(SystemLength/2)
                        if (min(c1,c2) <= SystemLength/2 and max(c1,c2) > SystemLength/2) \
                        and abs(c1-c2) <= SystemLength/2:
                            #print "updating sector!"
                            sector = (sector + 1)%2
                        #print "Error state after fixing:" + str(chain)'''
                        
                        
                            
        
        
        
        
        
        NumParallel = NumOfLambdas / 2
        for NumPar in xrange(NumParallel):
            AboutToSwap = CorrectionSwaps[(IndexInCycle + NumPar*2*BaseCycleLength)%NumberOfSwaps]
            
            #This if statement checks to see if the protocol is about to swap a defect off one of the
            #measurement rails.  If it is, it doesn't do the swap (i.e., defects should accumulate on
            #the measurement rails)
            if (not ((AboutToSwap in Attractors and chain[AboutToSwap] == 1) or \
                  (AboutToSwap+1)%SystemLength in Attractors and chain[(AboutToSwap+1)%SystemLength] == 1)):
                chain = CSwap(chain, AboutToSwap)  
                

                

        sector = (sector + CheckSector(IsingChain,p,pl,pr))%2 
        
        PHopInter, PAnnInter, PCreInter = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
        ##print PHopInter, PAnnInter, PCreInter
        if (PHop != PHopInter or PAnn != PAnnInter or PCre != PCreInter):
            ProbabilityHasChanged = True
        IndexInCycle = (IndexInCycle+1)%NumberOfSwaps
        CycleTime+=CorrectionPeriod/NumberOfSwaps
        
    
    NoHops = (PHop == 0)
    #####print "At end of correction", chain
    ##print "Starttime: ",oldtime,"Candidate endtime:",newtime,"Cycle endtime:",oldtime+CycleTime
    ##print "New rate equations: ", ProbabilityHasChanged, "Nohops: ", NoHops, "Proceeded?", Proceed
    
    return chain, CycleTime, ProbabilityHasChanged, NoHops, Proceed, sector
    
    
def CorrectionProtocol(chain, oldtime, newtime, CorrectionRate, PHop, PAnn, PCre, CorrectionSwaps, sector):
    
    ##print "What"
    CycleTime = 0
    #PHop, PAnn, PCre = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
    ##print PHop, PAnn, PCre
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
        #print oldtime + CycleTime
        #print "Chainstate: " + str(chain)
        p = int(floor(len(chain)/2.))
        pl,pr = chain[p],chain[p+1] #previous values of chain
        #####print "Timing information: ", CycleTime,"\t",oldtime,"\t", newtime,"\t",(newtime-oldtime)-CycleTime,"\t",CorrectionPeriod/NumberOfSwaps
        #chain = CSwap(chain, CorrectionSwaps[IndexInCycle])        
        #parallel?
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + NumberOfSwaps/4)%NumberOfSwaps])
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 2*NumberOfSwaps/4)%NumberOfSwaps])
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 3*NumberOfSwaps/4)%NumberOfSwaps])
        
        #NumParallel = 0       
        #if (SystemLength/3) % 2 == 0:
        #    NumParallel = SystemLength / 6
        #else:
        #NumParallel = int(floor(SystemLength/6.))
        NumParallel = NumOfLambdas / 2
        for NumPar in xrange(NumParallel):
            AboutToSwap = CorrectionSwaps[(IndexInCycle + NumPar*2*BaseCycleLength)%NumberOfSwaps]
            
            chain = CSwap(chain, AboutToSwap)  


 
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 2*NumberOfSwaps/4)%NumberOfSwaps])
        #chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 3*NumberOfSwaps/4)%NumberOfSwaps])
 

        sector = (sector + CheckSector(IsingChain,p,pl,pr))%2 
        
        PHopInter, PAnnInter, PCreInter = CalculateProbabilities(chain, CreationRate, AnnihilationRate, HoppingRate)
        ##print PHopInter, PAnnInter, PCreInter
        if (PHop != PHopInter or PAnn != PAnnInter or PCre != PCreInter):
            ProbabilityHasChanged = True
        IndexInCycle = (IndexInCycle+1)%NumberOfSwaps
        CycleTime+=CorrectionPeriod/NumberOfSwaps
        
    
    NoHops = (PHop == 0)
    #####print "At end of correction", chain
    ##print "Starttime: ",oldtime,"Candidate endtime:",newtime,"Cycle endtime:",oldtime+CycleTime
    ##print "New rate equations: ", ProbabilityHasChanged, "Nohops: ", NoHops, "Proceeded?", Proceed
    
    return chain, CycleTime, ProbabilityHasChanged, NoHops, Proceed, sector

def CheckState(chain, sector):
    if sum(chain)==0:
        return 2*sector-1
    else:
        return 0


def ProcessTraj(traj,avgtraj,maxtime):
    ##print traj
    #avgtraj[0]+=traj[0][1]
    trajindex = 0
    for i, val in enumerate(avgtraj):
        
        if i < len(traj) and i > 0: #safety first!
            
            while trajindex < len(traj) and traj[trajindex][0] < (1.0*maxtime / len(avgtraj))*i:
                ##print "Window: ",(1.0*maxtime / len(avgtraj))*i
                trajindex+=1
            avgtraj[i]+=traj[trajindex-1][1]
    
    return avgtraj






#Tests = SwapProtocol2(SystemLength)
##print Tests



#CorrectionRate = AnnihilationRate/len(SwapProtocol(SystemLength))/10.
#CorrectionRate = .001


#****
#CorrectionRate = .15*(1.5)**(-(float(sys.argv[1])-1))
#CorrectionRate = .0012
#CorrectionRate = .0015 * (1.094468)**(-(float(arg1)-1))
#CRE = .038085*(1. / SystemLength) - .00002055 #Empirically determined maximum correction rate
#CRE = 0.0000005156*(NumOfLambdas)**2.0 - 0.0000085000*NumOfLambdas + 0.0000940000


########################################################################Old params which work for L=96
#CRE = 0.0000002432*(NumOfLambdas)**2.0 + 0.0000029014*NumOfLambdas - 0.0000165727
#CorrectionRate = 4.*CRE * (1.094468)**(-(float(arg1)-1))
#CorrectionRate = 4.*CRE * (1.094468)**(-(float(15)-1))

CorrectionRate = float(sys.argv[2])

##print CRE
#print CorrectionRate
#print CreationRate
#This should put the maximum rate somewhere around the middle of the 30 jobs (because the 1.09 etc. factor is approx 1/4 when arg1 = 15)


#****

#CorrectionSwaps = SwapProtocol5(SystemLength)
#CorrectionSwaps = SwapProtocol6(NumOfLambdas)



##print SwapProtocol3(SystemLength)
##print CorrectionSwaps
zeros = []
ones = []
trajs = [0 for i in xrange(100)]
states = []
fliptimes = []
seed = numpy.random.rand()

NumOfLambdasList = [4,6,8,10,12,14] 

##############################################################Start of loop
for iii in xrange(10):
	##print i
	#IsingChain[0] = 1
	#IsingChain[5] = 1
	NumOfLambdas = NumOfLambdasList[int(sys.argv[1])%len(NumOfLambdasList)]
	#NumOfLambdas = 3
	SystemLength = (3+2*2)*NumOfLambdas
	IsingChain = [0 for i in xrange(SystemLength)]
	DefectAge = [0 for k in xrange(NumOfLambdas)]

	CorrectionSwaps = SwapProtocol7(SystemLength, SystemLength/NumOfLambdas)
	BaseCycleLength = len(CorrectionSwaps)/NumOfLambdas

	#Temperature = .12
	Temperature = .08
	IsingChain = [0 for i in xrange(SystemLength)]
	Delta = 1.0
	CreationRate = abs(1./(1-exp(Delta*1.0/Temperature)))
	AnnihilationRate = abs(1./(1-exp(-Delta*1.0/Temperature)))
	HoppingRate = Temperature

	counter = 0
	sector = 0
	totalt = 0

	traj = []
	state = []
	#while(counter < 1000 and totalt < exp(1./Temperature)):
	FlipHasOccurred = False


	while(FlipHasOccurred == False):# and counter < 100000:# and totalt < exp(1./Temperature)):
		#print "Overcounter: " + str(counter)
		curstate = CheckState(IsingChain,sector)
		if curstate == 1:
			FlipHasOccurred = True
			fliptimes.append(totalt)
		#if totalt >= 4*9955416204 and curstate != 1:
		#    FlipHasOccurred = True
		#    fliptimes.append(totalt)
		#    #print "Doublelifetime, ending"


		p = int(floor(len(IsingChain)/2.))
		pl,pr = IsingChain[p],IsingChain[p+1] #previous values of chain
		##print "Before correction", IsingChain
		IsingChain, t, sector = AdvanceTime(IsingChain, totalt, CreationRate, AnnihilationRate, HoppingRate, sector)
		#print IsingChain

		totalt+=t

		##print sector
		##print totalt, sector
		'''
		traj.append([totalt, sector])
		if sector == 1:
			ones.append(totalt)
		else:
			zeros.append(totalt)

		state.append([t,curstate])
		'''    
		counter+=1
		#print "Flip? " + str(FlipHasOccurred) + " time: " + str(totalt)
		#print "Flipchain" + str(IsingChain)
	
	with open('/home/daniel.freeman/SparseMeasurements/Runs/' + str(seed) + '.' + str(Temperature) + '.' + str(iii) + '.sweep.dat','a+') as f:
		f.write(str(Temperature)+"\t"+str(NumOfLambdas)+"\t"+str(SystemLength)+"\t"+str(CorrectionRate)+"\t"+str(totalt)+"\t"+str(totalt*totalt)+"\n")
#trajs = ProcessTraj(traj, trajs, 40000.)
#states.append(state)
#if FlipHasOccurred == False:
#	overcounts+=1
#	fliptimes.append(totalt)
###############################################End of loop


##print fliptimes
##print states[0]
#timeweighted = [[s[0]*s[1] for s in ss] for ss in states]
#avgstates = [float(sum(s))/len(s) for s in timeweighted]

##print trajs
##print avgstates
#a,b=histogram(ones, linspace(100,20000.,21))
##print "done"

print Temperature,"\t",NumOfLambdas,"\t",SystemLength,"\t",CorrectionRate,"\t",numpy.mean(fliptimes),"\t",numpy.std(fliptimes)/numpy.sqrt(len(fliptimes))
