import sys
import argparse as ap
import csv
import random
from collections import deque
from ContinuousTime import *
import resource

#resource.setrlimit(resource.RLIMIT_VMEM, (megs * 6000,megs * 8000))
resource.setrlimit(resource.RLIMIT_AS, (4500000000, 4500000100))

NontrivialRecombinationAverage = []
TimeToRecombineAverage = []

def PerformRun(Dim, Temp, Delt, TranRate, CurTime, TotTime, NumOfRuns, NumOfPoints, Mod, FilNam):


    InitialListOfExcitations = []
    
    for i in range(Dim):
        SubList = []
        for j in range(Dim):
            SubList.append(0)
        InitialListOfExcitations.append(SubList)
    
    
    Averages = []
    Deviations = []
    
    
    Delta = Delt
    Temperature = Temp
    SpawnRate = abs(1./(1-exp(Delta*1.0/Temperature)))
    AnnihilateRate = abs(1./(1-exp(-Delta*1.0/Temperature)))
    TranslateRate = Temperature
    CurrentTime = CurTime
    TotalTime = TotTime
    #TotalTime = 100
    NumberOfRuns = NumOfRuns
    Mode = Mod #can be Fine or Coarse.  Fine will use every timepoint in the plot, Coarse will look at discretized timepoints, effectively binning.
    NumberOfPoints = NumOfPoints
    CoarsePoints = linspace(0,TotalTime,NumberOfPoints)
    MyLattice = Lattice(InitialListOfExcitations,[0,0])
    
    TimesList = []
    SetOfTimesAndObservables = []
    DataSet = [] #stores [number of runs, [order params],[histogram],[windingslist]]
    
    def AccumulateObservables(Obs,SetOfTimesAndObservables):
        '''
	totalaveragelifetime = 0
        totalescapedlifetime = 0
        totalorderparam = 0.
        escapedlifetimecounter = 0
        allcounter = 0
        ecounter = 0
        totalsingleslifetime = 0.
        scounter = 0
        totalescapedsteps = 0
        estepscounter = 0
        totaldistance = 0.
        
        nontriv = 0
        averagelifetime = 0
        escapedlifetime = 0
        
        OrderParams = []
        
        
        def ExtractDistance(pair):
            return math.sqrt((pair[0][0]-pair[1][0])**2 + (pair[0][1]-pair[1][1])**2)
        
        
        #for distanceindex in xrange(len(SetOfTimesAndObservables[n][4])):
            
        
        for lifetimeindex in xrange(len(SetOfTimesAndObservables[3])/2-1):
            if round(SetOfTimesAndObservables[3][2*lifetimeindex][0],5)!=round(SetOfTimesAndObservables[3][2*lifetimeindex+1][0],5):
                nontriv+=1
        
        for lifetimeindex in xrange(len(SetOfTimesAndObservables[3])):
            allcounter+=1
            totalaveragelifetime+=SetOfTimesAndObservables[3][lifetimeindex][0]
            
            if SetOfTimesAndObservables[3][lifetimeindex][1] != 1:
                ecounter+=1
                totalescapedlifetime+=float(SetOfTimesAndObservables[3][lifetimeindex][0])
                estepscounter+=1
                totalescapedsteps+=float(SetOfTimesAndObservables[3][lifetimeindex][1])
                
                totaldistance+=ExtractDistance(SetOfTimesAndObservables[4][lifetimeindex])
                
            else:
                totalsingleslifetime+=float(SetOfTimesAndObservables[3][lifetimeindex][0])
                scounter+=1
        
        #averagelifetime = averagelifetime / (1.0*len(SetOfTimesAndObservables[n][3]))
        #totalaveragelifetime += averagelifetime
        totalorderparam += nontriv/(len(SetOfTimesAndObservables[3])/2.)
        #totalescapedlifetime += escapedlifetime
        #print averagelifetime
        
        #print "Average lifetime:\t" + str(totalaveragelifetime / (1.0*allcounter))
        #print "Average order parameter:\t" + str(totalorderparam / (1.0*NumberOfRuns))
        #print "Average escaped lifetime:\t" + str(totalescapedlifetime / (1.0*ecounter))
        #print "Average singles lifetime:\t" + str(totalsingleslifetime / (1.0*scounter))
        #print "Average number of escaped steps:\t" + str(totalescapedsteps / (1.0*estepscounter))
        #print "Average escaped distance:\t" + str(totaldistance / (1.0*estepscounter))
        #print "All:\t" + str(allcounter)
        #print "Esc:\t" + str(ecounter)
        
        
        
        
        #print [n[3] for n in SetOfTimesAndObservables]
        partition = numpy.linspace(1,1000,1000)
        
        stepdata = [float(q[1]) for q in SetOfTimesAndObservables[3]]
        #print partition
        #print [entry[1] for q in [n[3] for n in SetOfTimesAndObservables] for entry in q]
        thishist = numpy.histogram(stepdata,partition)
        #print "Number of walks:\t" + str(len(stepdata))
        #print thishist[0]
        #os.chdir("C:\Users\Daniel\Dropbox\EclipseWorkspace\MyToricCodeProject\src\Histograms")
        
        #numpy.savetxt(FilNam[:-3]+".hist.csv",thishist[0],delimiter=",")
        if allcounter == 0:
            allcounter = 1
        if ecounter == 0:
            ecounter = 1
        if scounter == 0:
            scounter = 1
        if estepscounter == 0:
            estepscounter = 1        
        
        OrderParams = [totalaveragelifetime / (1.0*allcounter),totalorderparam / (1.0*NumberOfRuns),totalescapedlifetime / (1.0*ecounter),\
            totalsingleslifetime / (1.0*scounter),totalescapedsteps / (1.0*estepscounter),totaldistance / (1.0*estepscounter)\
            ,allcounter,ecounter,len(stepdata)]
        
	#*******************************************************************************************************************
        '''
        
	OrderParams = SetOfTimesAndObservables[3]

        
        
        
        CurrentIndexes = 0
        counter = 0
        ObservablesList = []

        for FineTime in CoarsePoints:
            
            #print "Writing " + str(counter) + " out of " + str(len(CoarsePoints))
            counter += 1
            Winding00 = 0
            Winding01 = 0
            Winding10 = 0
            Winding11 = 0
            
            NumOfExcitations = 0
            
            for n in xrange(NumberOfRuns):
                #print "Current indexes [n]: " + str(CurrentIndexes)
                #NumOfExcitations=0
                NumOfExcitations+=SetOfTimesAndObservables[2][CurrentIndexes]
                
                '''if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,0] or \
                SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,1]:
                    WindingX1 +=1
                if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,0] or \
                SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,1]:
                    WindingX2 +=1'''
                
                if SetOfTimesAndObservables[1][CurrentIndexes] == [0,0]:
                    Winding00 +=1
                if SetOfTimesAndObservables[1][CurrentIndexes] == [1,0]:
                    Winding10 +=1
                if SetOfTimesAndObservables[1][CurrentIndexes] == [0,1]:
                    Winding01 +=1
                if SetOfTimesAndObservables[1][CurrentIndexes] == [1,1]:
                    Winding11 +=1
                
            
            NumOfExcitations = NumOfExcitations*1.0 / NumberOfRuns*1.0
            
            Winding00 = Winding00*1.0 / NumberOfRuns*1.0
            Winding10 = Winding10*1.0 / NumberOfRuns*1.0
            Winding01 = Winding01*1.0 / NumberOfRuns*1.0
            Winding11 = Winding11*1.0 / NumberOfRuns*1.0

            
            ObservablesList.append([FineTime,Winding00,Winding01,Winding10,Winding11,NumOfExcitations,NumOfExcitations**2])
            
            
            #This inner loop only goes to len(this specific run) - 1 because it always checks the next
            #entry to see if it's a time after FineTime
            CoarseTime = CurrentIndexes
            while CoarseTime < len(SetOfTimesAndObservables[0]) \
            and SetOfTimesAndObservables[0][CoarseTime] <= FineTime:
                CoarseTime +=1
            CurrentIndexes = CoarseTime - 1            
        #CurrentIndexes is properly keeping track of the correct time to check (by its index) in each of the
        #independent runs.
            
        #print SetOfTimesAndObservables[0][2]
        if len(Obs)==0:
            Obs.append(0)
            #Obs.append([totalaveragelifetime / (1.0*allcounter),totalorderparam / (1.0*NumberOfRuns),totalescapedlifetime / (1.0*ecounter),\
            #totalsingleslifetime / (1.0*scounter),totalescapedsteps / (1.0*estepscounter),totaldistance / (1.0*estepscounter)\
            #,allcounter,ecounter,len(stepdata)])
            #Obs.append(list(thishist[0]))
            Obs.append(OrderParams)
	    Obs.append(ObservablesList)
            #print Obs[3]
            
        else:
            CurrentIndex = Obs[0]
            #print Obs[3][-10:]
            #print ObservablesList[-10:]
            Obs[0] = Obs[0] + 1
            Obs[1] = [(n + o*CurrentIndex)/(CurrentIndex+1.) for n,o in zip(OrderParams,Obs[1])]
            #Obs[2] = [(n + o*CurrentIndex)/(CurrentIndex+1.) for n,o in zip(thishist[0],Obs[2])]
            Obs[2] = [[(n + o*CurrentIndex)/(CurrentIndex+1.) for n,o in zip(ObservablesList[k],Obs[2][k])] for k in xrange(len(Obs[2]))]
        
        return Obs
    
    

    
    for RunNumber in xrange(NumberOfRuns):
        print "Run number: " + str(RunNumber)
        MyLattice = Lattice(InitialListOfExcitations,[0,0])
        #Times, Observables, Occupations, Lifetimes, Distances = Simulate3(MyLattice, SpawnRate, AnnihilateRate, TranslateRate, CurrentTime, TotalTime)
        Times, Observables, Occupations, OrderPs = Simulate3(MyLattice, SpawnRate, AnnihilateRate, TranslateRate, CurrentTime, TotalTime)

	DataSet = AccumulateObservables(DataSet,[Times,Observables,Occupations,OrderPs])

	
	#SetOfTimesAndObservables.append([Times,Observables,Occupations,Lifetimes,Distances])
        

	#DataSet = AccumulateObservables(DataSet,[Times,Observables,Occupations,Lifetimes,Distances])
        #print Lifetimes
        #print sum(sorted(Lifetimes, key = lambda x: x[0]))
        for t in xrange(len(Times)):
            TimesList.append([t,Times[t],RunNumber])
       
        if RunNumber%10 == 0:
            with open(FilNam,"w+") as f:
                f.truncate()    
                f.write("\t".join(str(o) for o in DataSet[1])+"\t"+str(DataSet[0])+"\n")
                for O in DataSet[2]:
                    f.write("\t".join([str(p) for p in O])+"\n")
            '''
	    with open(FilNam[:-3]+"hist.dat","w+") as f:
                f.truncate()    
                for O in DataSet[2]:
                    f.write(str(O)+"\n")
	    '''
    
    
    
    
    
    #This loop picks an element of t.  It then finds the time points in each of the other
    #runs that are closest to t.  It keeps track of where it is, relatively, in each of the 
    #sublists for the times of each independent runs with CurrentIndex list
    DataTimes = []
    DataAverages = []
    
    
    return DataSet
    
    
tempdimtimes = []

#with open('PitzerLargeParameterSet.dat','r') as f:
#with open('BetweenParameterSet.dat','r') as f:
#with open('HighTFits.out','r') as f:
with open('BetweenParameterSet2.dat','r') as f:
	for line in f:
		coeffs = line.split("\t")
		tempdimtimes.append([float(coeffs[0]),int(coeffs[1]),float(coeffs[4])])
		                    
#for LargeRun in xrange(3):
for LargeRun in xrange(1):   
    #Times = [250000,25000,5000,2500,1800,1500,1100,700,300]
    #TArray = [.07,.08,.09,.1,.11,.12,.13,.14,.15]
    
    #print "Weee"
    #Times = [250000,25000,5000,2500,1800,1500,1100,700,300]
    #TArray = [.07,.08,.09,.1,.11,.12,.13,.14,.15]
    
    #print "Weee"
    
    Dim = 8*(2**(LargeRun+1))
    #Dim = 4

    Temp = .06
    #Temp = TArray[LargeRun]
    
    
    Delt = 1.
    TranRate = Temp
    CurTime = 0
    
    #TotTime = 100
    #TotTime = Times[LargeRun]
    TotTime = 22000
    
    #NumOfRuns = 200
    NumOfRuns = 1000
    
    NumOfPoints = 1000
    SubRange = 3
    Mod = "Coarse"
    #TArray = [.08,.09,.1]
    #Times = [[11500,3140,1200],[4000,1100,440],[1900,670,343]]
    TArray = [.11,.12,.13,.14,.15]
    Times = [[430,220,150],[200,110,80],[100,60,50],[60,40,30],[35,25,20]]    

    #counter = 0
    #for SubRun in xrange(SubRange):
    parameterindex = (int(sys.argv[3])-1)*6
    #parameterindex = int(sys.argv[3])%20
    for SubRun in xrange(10):
        #TotTime = ((SubRun+1)*100)*TotTime/(10.*(LargeRun+1.))
        #Temp = .07 - SubRun*.01
        #TotTime = Times[LargeRun][SubRun]
    	#Temp = TArray[SubRun]
    	
        #SubRun = 542+counter
        #counter+=1
	SubRunI = (parameterindex + SubRun)%len(tempdimtimes)

        Dim = tempdimtimes[SubRunI][1]
    	Temp = tempdimtimes[SubRunI][0]
    	TotTime = tempdimtimes[SubRunI][2]
        TranRate = Temp
        #print Dim
        #print Temp
        #print "Subrun number: " + str(SubRun)
        FilNam = "MeasuredRuns9/" + str(Dim) + "_Temp" + str(Temp) + "_Delt" + str(Delt) + "_TransRate" + str(TranRate) + "_" +\
        str(LargeRun) + "_" + str(SubRun) + "_" + str(sys.argv[2]) + "_" + str(sys.argv[3])+ ".csv"
        Observables = []
        #for x in xrange(NumOfRuns):
        #print x
        Observables = PerformRun(Dim, Temp, Delt, TranRate, CurTime, TotTime, NumOfRuns, NumOfPoints, Mod, FilNam)
        
        
        #print "\t".join(str(o) for o in Observables[1])
        
        
    
    #for i in xrange(LargeRun*SubRange, (LargeRun+1)*SubRange):
        #print str(NontrivialRecombinationAverage[i][0]) + "," + str(NontrivialRecombinationAverage[i][1]) + "," + str(NontrivialRecombinationAverage[i][2])
        #print str(TimeToRecombineAverage[i][0]) + "," + str(TimeToRecombineAverage[i][1]) + "," + str(TimeToRecombineAverage[i][2])
    
    
    
    
    
    '''times = []
    
    for i in xrange(2000):
        times.append(0)
        
    for SubRun in xrange(SubRange):
        FilNam = "src/runs/Summer2013Runs/TrivialRecombination/" + str(Dim) + "_Temp" + str(Temp) + "_Delt" + str(Delt) + "_TransRate" + str(TranRate) + "_" +\
        str(LargeRun) + "_" + str(SubRun) + ".csv"
        
        with open(FilNam, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            index = 0
            for row in reader:
                
                if SubRun == 0:
                    times[index]=float(row[0])
                
                print row[1]
                windings[index]+=float(row[1])-.5
                #deviations[index]+=(float(row[1])-.5)**2
                index+=1
    
    for w in windings:
        print str(w/89.)'''
        
    
    
    '''deviations = []
    times = []
    windings = []
    for i in xrange(2000):
        times.append(0)
        windings.append(0)
        deviations.append(0)
    
    for SubRun in xrange(SubRange):
        
        FilNam = "runs/testbatch/Dim" + str(Dim) + "_Temp" + str(Temp) + "_Delt" + str(Delt) + "_TransRate" + str(TranRate) + "_" +\
        str(SubRun) + ".csv"
        
        with open(FilNam, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            index = 0
            for row in reader:
                
                if SubRun == 0:
                    times[index]=float(row[0])
                
                print row[1]
                windings[index]+=float(row[1])-.5
                #deviations[index]+=(float(row[1])-.5)**2
                index+=1
    
    for w in windings:
        print str(w/89.)'''


            
        
        
        
        
        
        
