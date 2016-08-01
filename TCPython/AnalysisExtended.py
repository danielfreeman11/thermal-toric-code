import sys
import argparse as ap
import csv
import random
from collections import deque
from ContinuousTime import *
import resource

#resource.setrlimit(resource.RLIMIT_VMEM, (megs * 6000,megs * 8000))
resource.setrlimit(resource.RLIMIT_AS, (3000000000, 3000000100))

NontrivialRecombinationAverage = []
TimeToRecombineAverage = []

def PerformRun(Dim, Temp, Delt, TranRate, CurTime, TotTime, NumOfRuns, NumOfPoints, Mod, FilNam):


    InitialListOfExcitations = []
    
    for i in range(Dim):
        SubList = []
        for j in range(Dim):
            SubList.append(0)
        InitialListOfExcitations.append(SubList)
    
    #sample = random.randint(0,2)
    #choices = [(2,4),(4,4),(3,5)]
    #choicex,choicey = choices[sample]
    ###########
    #InitialListOfExcitations[3][3] = 1
    #InitialListOfExcitations[choicex][choicey] = 1
    ###########
    
    
    Averages = []
    Deviations = []
    
    
    Delta = Delt
    Temperature = Temp
    #SpawnRate = exp(-Delta/2./Temperature)
    #SpawnRate = 0
    #AnnihilateRate = exp(Delta/2./Temperature)
    #TranslateRate = TranRate
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
    
    
    
    #Delta = 1.
    #Temperature = 1/14.
    #SpawnRate = exp(-Delta/2/Temperature)
    #AnnihilateRate = exp(Delta/2/Temperature)
    #TranslateRate = 2
    #CurrentTime = 0
    #TotalTime = 10000
    #TotalTime = 100
    #NumberOfRuns = 1000
    #Mode = "Coarse" #can be Fine or Coarse.  Fine will use every timepoint in the plot, Coarse will look at discretized timepoints, effectively binning.
    #NumberOfPoints = 1000
    #CoarsePoints = linspace(0,TotalTime,NumberOfPoints)
    
    MyLattice = Lattice(InitialListOfExcitations,[0,0])
    
    #Times, Observables = Simulate2(MyLattice, SpawnRate, AnnihilateRate, TranslateRate, CurrentTime, TotalTime)
    
    #for t in xrange(len(Times)): print str(Times[t]) + "\t" + str(Observables[t])
    
    
    
    
    TimesList = []
    SetOfTimesAndObservables = []
    
    

    
    for RunNumber in xrange(NumberOfRuns):
        #print "Run number: " + str(RunNumber)
        MyLattice = Lattice(InitialListOfExcitations,[0,0])
        Times, Observables, Occupations, Lifetimes, Distances = Simulate3(MyLattice, SpawnRate, AnnihilateRate, TranslateRate, CurrentTime, TotalTime)
        #Times, Observables, Occupations = SimulateSingleShot(MyLattice, 0.0, AnnihilateRate, TranslateRate, CurrentTime, TotalTime)
        SetOfTimesAndObservables.append([Times,Observables,Occupations,Lifetimes,Distances])
        #print Lifetimes
        #print sum(sorted(Lifetimes, key = lambda x: x[0]))
        for t in xrange(len(Times)):
            TimesList.append([t,Times[t],RunNumber])

    
    
    
    
    
    #This loop picks an element of t.  It then finds the time points in each of the other
    #runs that are closest to t.  It keeps track of where it is, relatively, in each of the 
    #sublists for the times of each independent runs with CurrentIndex list
    DataTimes = []
    DataAverages = []
    
    if(Mode == "TrivialSingleShot"):
        EndTimesAverage = 0
        for Run in xrange(NumberOfRuns):
            EndTimesAverage += SetOfTimesAndObservables[Run][0][-1]

                
        EndTimesAverage = EndTimesAverage / (NumberOfRuns*1.0)
        TimeToRecombineAverage.append((double(EndTimesAverage),double(Dim),double(NumberOfRuns)))
        print EndTimesAverage
    
    
    if(Mode == "SingleShot"):

        NontrivialRecombinations = 0
        for Run in xrange(NumberOfRuns):
            EndTime = SetOfTimesAndObservables[Run][0][-1]

            WindingX1 = 0
            WindingX2 = 0
            
            NumOfExcitations = 0
            
            n=Run
            NumOfExcitations+=SetOfTimesAndObservables[n][2][-1]
            
            if SetOfTimesAndObservables[n][1][-1] == [0,0] or \
            SetOfTimesAndObservables[n][1][-1] == [0,1]:
                WindingX1 +=1
            if SetOfTimesAndObservables[n][1][-1] == [1,0] or \
            SetOfTimesAndObservables[n][1][-1] == [1,1]:
                WindingX2 +=1
                NontrivialRecombinations += 1
                
        NontrivialRecombinations = NontrivialRecombinations / (NumberOfRuns*1.0)
        NontrivialRecombinationAverage.append((double(NontrivialRecombinations),double(Dim),double(NumberOfRuns)))
        #print (NontrivialRecombinations)

                
        '''with open(FilNam, "w") as outFile:
            writer = csv.writer(outFile)
            row = deque([WindingX1,WindingX2,NumOfExcitations])
            row.appendleft(double(EndTime))
            writer.writerow(row)
            outFile.flush()'''
        
        
    
    
    if(Mode == "Coarse"):
        #print("Performing a coarse filesave")
        
		
	#***************************************************************************************************
        #nontriv = 0
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
        
        for n in xrange(NumberOfRuns):
            nontriv = 0
            averagelifetime = 0
            escapedlifetime = 0
            
            
            def ExtractDistance(pair):
                return math.sqrt((pair[0][0]-pair[1][0])**2 + (pair[0][1]-pair[1][1])**2)
            
            
            #for distanceindex in xrange(len(SetOfTimesAndObservables[n][4])):
                
            
            for lifetimeindex in xrange(len(SetOfTimesAndObservables[n][3])/2-1):
                if round(SetOfTimesAndObservables[n][3][2*lifetimeindex][0],5)!=round(SetOfTimesAndObservables[n][3][2*lifetimeindex+1][0],5):
                    nontriv+=1
            
            for lifetimeindex in xrange(len(SetOfTimesAndObservables[n][3])):
                allcounter+=1
                totalaveragelifetime+=SetOfTimesAndObservables[n][3][lifetimeindex][0]
                
                if SetOfTimesAndObservables[n][3][lifetimeindex][1] != 1:
                    ecounter+=1
                    totalescapedlifetime+=float(SetOfTimesAndObservables[n][3][lifetimeindex][0])
                    estepscounter+=1
                    totalescapedsteps+=float(SetOfTimesAndObservables[n][3][lifetimeindex][1])
                    
                    totaldistance+=ExtractDistance(SetOfTimesAndObservables[n][4][lifetimeindex])
                    
                else:
                    totalsingleslifetime+=float(SetOfTimesAndObservables[n][3][lifetimeindex][0])
                    scounter+=1
            
            #averagelifetime = averagelifetime / (1.0*len(SetOfTimesAndObservables[n][3]))
            #totalaveragelifetime += averagelifetime
            totalorderparam += nontriv/(len(SetOfTimesAndObservables[n][3])/2.)
            #totalescapedlifetime += escapedlifetime
            #print averagelifetime
        
        print "Average lifetime:\t" + str(totalaveragelifetime / (1.0*allcounter))
        print "Average order parameter:\t" + str(totalorderparam / (1.0*NumberOfRuns))
        print "Average escaped lifetime:\t" + str(totalescapedlifetime / (1.0*ecounter))
        print "Average singles lifetime:\t" + str(totalsingleslifetime / (1.0*scounter))
        print "Average number of escaped steps:\t" + str(totalescapedsteps / (1.0*estepscounter))
        print "Average escaped distance:\t" + str(totaldistance / (1.0*estepscounter))
        print "All:\t" + str(allcounter)
        print "Esc:\t" + str(ecounter)
        
        
        #print [n[3] for n in SetOfTimesAndObservables]
        partition = numpy.linspace(1,1000,1000)
        stepdata = [float(entry[1]) for q in [n[3] for n in SetOfTimesAndObservables] for entry in q]
        #print partition
        #print [entry[1] for q in [n[3] for n in SetOfTimesAndObservables] for entry in q]
        thishist = numpy.histogram(stepdata,partition)
        print "Number of walks:\t" + str(len(stepdata))
        #print thishist[0]
        #os.chdir("C:\Users\Daniel\Dropbox\EclipseWorkspace\MyToricCodeProject\src\Histograms")
        
        numpy.savetxt(FilNam[:-3]+".hist.csv",thishist[0],delimiter=",")
        
        
        #*******************************************************************************************************************
        
        CurrentIndexes = []
        for i in xrange(NumberOfRuns):
            CurrentIndexes.append(0)
        
        counter = 0
        #with open("data/test" + str(time.time()) + ".csv", "w") as outFile:
        with open(FilNam, "w") as outFile:
            writer = csv.writer(outFile)
            for FineTime in CoarsePoints:
                
                #print "Writing " + str(counter) + " out of " + str(len(CoarsePoints))
                counter += 1
                Winding00 = 0
                Winding01 = 0
                Winding10 = 0
                Winding11 = 0
                
                NumOfExcitations = 0
                
                for n in xrange(NumberOfRuns):
                    #print "Current indexes [n]: " + str(CurrentIndexes[n])
                    #NumOfExcitations=0
                    NumOfExcitations+=SetOfTimesAndObservables[n][2][CurrentIndexes[n]]
                    
                    '''if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,0] or \
                    SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,1]:
                        WindingX1 +=1
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,0] or \
                    SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,1]:
                        WindingX2 +=1'''
                    
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,0]:
                        Winding00 +=1
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,0]:
                        Winding10 +=1
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,1]:
                        Winding01 +=1
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,1]:
                        Winding11 +=1
                    
                
                NumOfExcitations = NumOfExcitations*1.0 / NumberOfRuns*1.0
                
                Winding00 = Winding00*1.0 / NumberOfRuns*1.0
                Winding10 = Winding10*1.0 / NumberOfRuns*1.0
                Winding01 = Winding01*1.0 / NumberOfRuns*1.0
                Winding11 = Winding11*1.0 / NumberOfRuns*1.0
                
                
                '''WindingX1Deviation = 0
                WindingX2Deviation = 0
                
                for n in xrange(NumberOfRuns):
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,0] or \
                    SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,1]:
                        WindingX1Deviation += (WindingX1 - 1)**2
                    else:
                        WindingX1Deviation += (WindingX1 - 0)**2
                    
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,0] or \
                    SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,1]:
                        WindingX2Deviation += (WindingX2 - 1)**2
                    else:
                        WindingX2Deviation += (WindingX2 - 0)**2
                
                WindingX1Deviation = WindingX1Deviation / NumberOfRuns
                WindingX2Deviation = WindingX2Deviation / NumberOfRuns
                
                WindingX1Deviation = sqrt(WindingX1Deviation)
                WindingX2Deviation = sqrt(WindingX2Deviation)'''
                
                
                #DataTimes.append(FineTime[1])
                #DataAverages.append([AveragedObservable,FineTime[0]])
                
                row = deque([Winding00,Winding01,Winding10,Winding11,NumOfExcitations])
                row.appendleft(FineTime)
                writer.writerow(row)
                outFile.flush()
                
                
                for n in xrange(NumberOfRuns):
                    #This inner loop only goes to len(this specific run) - 1 because it always checks the next
                    #entry to see if it's a time after FineTime
                    CoarseTime = CurrentIndexes[n]
                    while CoarseTime < len(SetOfTimesAndObservables[n][0]) \
                    and SetOfTimesAndObservables[n][0][CoarseTime] <= FineTime:
                        '''print "***"
                        print "CoarseTime: " + str(CoarseTime)
                        print "Actual Time: " + str(SetOfTimesAndObservables[n][0][CoarseTime])
                        print "Finetime: " + str(FineTime)
                        print "***"'''
                        CoarseTime +=1
                    CurrentIndexes[n] = CoarseTime - 1            
            #CurrentIndexes is properly keeping track of the correct time to check (by its index) in each of the
            #independent runs.
        
    #print SetOfTimesAndObservables[0][2]
    
    
    
    
    if(Mode == "Fine"):
        TimesList = sorted(TimesList, key = lambda x: x[1])
        # (TimesList)
        
        CurrentIndexes = []
        for i in xrange(NumberOfRuns):
            CurrentIndexes.append(0)
        
        counter = 0
        with open(FilNam, "w") as outFile:
            writer = csv.writer(outFile)
            for FineTime in TimesList:
                
                #print "Writing " + str(counter) + " out of " + str(len(TimesList))
                counter += 1
                WindingX1 = 0
                WindingX2 = 0
                
                NumOfExcitations = 0
                        
                for n in xrange(NumberOfRuns):
                    #print "Current indexes [n]: " + str(CurrentIndexes[n])
                    #NumOfExcitations=0
                    NumOfExcitations+=SetOfTimesAndObservables[n][2][CurrentIndexes[n]]
                    
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,0] or \
                    SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [0,1]:
                        WindingX1 +=1
                    if SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,0] or \
                    SetOfTimesAndObservables[n][1][CurrentIndexes[n]] == [1,1]:
                        WindingX2 +=1
                
                NumOfExcitations = NumOfExcitations*1.0 / NumberOfRuns*1.0
                
                WindingX1 = WindingX1*1.0 / NumberOfRuns*1.0
                WindingX2 = WindingX2*1.0 / NumberOfRuns*1.0
                
                
                #AveragedObservable = AveragedObservable*1.0 / NumberOfRuns*1.0
                #DataTimes.append(FineTime[1])
                #DataAverages.append([AveragedObservable,FineTime[0]])
                
                row = deque([WindingX1,WindingX2,NumOfExcitations])
                row.appendleft(FineTime[1])
                writer.writerow(row)
                outFile.flush()
                
                
                for n in xrange(NumberOfRuns):
                    #This inner loop only goes to len(this specific run) - 1 because it always checks the next
                    #entry to see if it's a time after FineTime
                    CoarseTime = CurrentIndexes[n]
                    while CoarseTime < len(SetOfTimesAndObservables[n][0]) \
                    and SetOfTimesAndObservables[n][0][CoarseTime] <= SetOfTimesAndObservables[FineTime[2]][0][FineTime[0]]:
                        CoarseTime +=1
                    CurrentIndexes[n] = CoarseTime - 1
            
            #CurrentIndexes is properly keeping track of the correct time to check (by its index) in each of the
            #independent runs.
            
            
'''for LargeRun in xrange(10):
    Dim = 10
    Temp = .07 - (LargeRun/(20.*(100./3.)))
    Delt = 1.
    TranRate = 1.
    CurTime = 0
    TotTime = 6000
    NumOfRuns = 1000
    NumOfPoints = 1000
    SubRange = 10
    Mod  = "Coarse"
    
    for SubRun in xrange(SubRange):
        
        FilNam = "runs/temperaturebatch6/Dim" + str(Dim) + "_Temp" + str(Temp) + "_Delt" + str(Delt) + "_TransRate" + str(TranRate) + "_" + str(SubRun)\
        + ".csv" 
        PerformRun(Dim, Temp, Delt, TranRate, CurTime, TotTime, NumOfRuns, NumOfPoints, Mod, FilNam)
        
    AnalysisFile = [[],[],[]]
    times = []
    windings = []
    
    deviations = []
    for i in xrange(1000):
        times.append(0)
        windings.append(0)
        deviations.append(0)
    
    for SubRun in xrange(SubRange):
        FilNam = "runs/temperaturebatch6/Dim" + str(Dim) + "_Temp" + str(Temp) + "_Delt" + str(Delt) + "_TransRate" + str(TranRate) + "_" + str(SubRun)\
        + ".csv"
    
        with open(FilNam, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            index = 0
            for row in reader:
                times[index]=float(row[0])
                windings[index]+=float(row[1])-.5
                deviations[index]+=(float(row[1])-.5)**2
                index+=1
                
    for i in xrange(len(deviations)):
        deviations[i] = math.sqrt(deviations[i]/float(SubRange) - (windings[i]/float(SubRange))**2)
        
    for i in xrange(len(windings)):
        windings[i] = windings[i]/float(SubRange)
    
    logwindings = []
    fittimes = []
    lastindex = 0
    while windings[lastindex] > 0 and lastindex < len(windings)-1:
        logwindings.append(math.log(windings[lastindex]))
        fittimes.append(times[lastindex])
        lastindex+=1
    
        
    #for i in xrange(len(deviations)):
     #   print str(times[i]) + "," + str(windings[i]) + "," + str(deviations[i])
                
    #print "Number of good points: " + str(lastindex)
    
    fittimes = numpy.array(fittimes)
    logwindings = numpy.array(logwindings)
    
    A = numpy.vstack([fittimes, numpy.ones(len(fittimes))]).T
    
    solution = numpy.linalg.lstsq(A, logwindings)
    m,c = solution[0]
    r = solution[1] 
    
    #print ("Slope: " + str(m) + " and offset: " + str(c) + " and residue: " + str(r))
    print(str(Temp) + "," + str(m) + "," + str(c))
    #print "what"
'''
tempdimtimes = []
with open('SweepParameters.dat','r') as f:
	for line in f:
		coeffs = line.split("\t")
		tempdimtimes.append([float(coeffs[0]),int(coeffs[1]),float(coeffs[4])])
		                    
#for LargeRun in xrange(3):
for LargeRun in xrange(50):   
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
    
    NumOfRuns = 10
    
    NumOfPoints = 1000
    SubRange = 3
    Mod = "Coarse"
    #TArray = [.08,.09,.1]
    #Times = [[11500,3140,1200],[4000,1100,440],[1900,670,343]]
    TArray = [.11,.12,.13,.14,.15]
    Times = [[430,220,150],[200,110,80],[100,60,50],[60,40,30],[35,25,20]]    


    #for SubRun in xrange(SubRange):
    for SubRun in xrange(len(tempdimtimes)):
        #TotTime = ((SubRun+1)*100)*TotTime/(10.*(LargeRun+1.))
        #Temp = .07 - SubRun*.01
        #TotTime = Times[LargeRun][SubRun]
	#Temp = TArray[SubRun]
	Dim = tempdimtimes[SubRun][1]
	Temp = tempdimtimes[SubRun][0]
	TotTime = tempdimtimes[SubRun][2]
	
	
	TranRate = Temp 
	print "Subrun number: " + str(SubRun)
        FilNam = "MeasuredRuns2/" + str(Dim) + "_Temp" + str(Temp) + "_Delt" + str(Delt) + "_TransRate" + str(TranRate) + "_" +\
        str(LargeRun) + "_" + str(SubRun)  + "_" + str(sys.argv[2]) + "_" + str(sys.argv[3])+ ".csv"
        PerformRun(Dim, Temp, Delt, TranRate, CurTime, TotTime, NumOfRuns, NumOfPoints, Mod, FilNam)
        
    
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


            
        
        
        
        
        
        
