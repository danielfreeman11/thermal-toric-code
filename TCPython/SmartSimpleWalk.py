from scipy import random
import math

Directions = []

Directions.append((0,1))
Directions.append((0,-1))
Directions.append((1,0))
Directions.append((-1,0))




#Walker will start at (0,2), (1,1), or ((-1%Dimx),1).  If it becomes adjacent to (0,0), it will annihilate.
#I've set the winding axis to be a couple columns away from the starting point for simplicity.
def SingleRun(Dimx,Dimy):

    def mapped(x,y):
        return ((x%Dimx),(y%Dimy))

    #For the calculation we want to perform, the dynamics look like a pair of random walkers appearing on a lattice in cells adjacent to one another.
    #The dynamics of these walkers allow them to randomly wiggle around the lattice, or annihilate each other (if they happen to be next to each other,
    #and if they happen to decide to annihilate--they don't necessarily annihilate if they are next to each other (although it is incredibly likely
    #they will annihilate one another if they're next to each other)).
    #This set of starting positions mimics the situation of a pair of walkers appearing and *NOT* immediately annihilating.  Thus, one of the walkers
    #moves either left, right, or up (you can see this by plotting the three starting positions).  If your lattice size was 10x10, (-1%Dimx) would
    #yield 9.
    StartingPositions = []
    StartingPositions.append(mapped(0,2))
    StartingPositions.append(mapped(0,-2))
    StartingPositions.append(mapped(2,0))
    StartingPositions.append(mapped(-2,0))
    StartingPositions.append(mapped(1,1))
    StartingPositions.append(mapped(-1,-1))
    StartingPositions.append(mapped(-1,1))
    StartingPositions.append(mapped(1,-1))

    
    StartingPositions.append(mapped(1,1))
    StartingPositions.append(mapped(-1,-1))
    StartingPositions.append(mapped(-1,1))
    StartingPositions.append(mapped(1,-1))
    
    #StartingPositions.append(((-1%Dimx),1))


    #We want to randomly sample from these three starting positions as they are all equally likely
    xPos,yPos = StartingPositions[random.randint(0,11)]
    
    #These are just some variables I care about for the simulation.
    xDistance,yDistance = 0,0
    lifetime = 0
    SectorX = 0
    SectorY = 0
    WindingAxisX = math.floor(Dimx/2.)
    WindingAxisY = math.floor(Dimy/2.)
    
    
    #A simple output statement to check that the random walker code works.
    def PrintState():
        for i in xrange(Dimx):
            print "\n"
            for j in xrange(Dimy):
                if (i,j) == (xPos,yPos):
                    print "1",
                else: print "0",
        print "\n"
    
    #This code advances everything I care about.  It takes the current position of the walker, and some other variables, and iterates them as necessary
    def Advance(xPos,yPos,xDistance,yDistance,SectorX,SectorY,lifetime):
        MoveX,MoveY = Directions[random.randint(0,4)]
        
        xDistance+=MoveX
        yDistance+=MoveY
        lifetime+=1
        
        if xPos == WindingAxisX and MoveX == -1:
            SectorX = (SectorX + 1)%2
        if xPos == WindingAxisX - 1 and MoveX == 1:
            SectorX = (SectorX + 1)%2
        
        if yPos == WindingAxisY and MoveY == -1:
            SectorY = (SectorY + 1)%2
        if yPos == WindingAxisY - 1 and MoveY == 1:
            SectorY = (SectorY + 1)%2
        
        xPos = (xPos + MoveX)%Dimx
        yPos = (yPos + MoveY)%Dimy
        
        return xPos,yPos,xDistance,yDistance,SectorX,SectorY,lifetime
    
    #This is the end condition for the random walk.  If the random walker comes within 1 unit of (0,0), the program ends because
    #it's overwhelming more likely for the random walkers to self-annihilate than they are to keep random walking.  Note that only one of the random
    #walkers is actually randomly walking in this simulation--one of them is pinned at 0,0.  This is valid because 2 random walkers walking on a lattice
    #can be thought of as just one random walker walking around, where we fix the origin on the other random walker.
    def CheckEnd(xPos,yPos):
        if ((xPos,yPos) == (1,0) or (xPos,yPos) == (0,1) or (xPos,yPos) == (Dimx-1,0) or (xPos,yPos) == (0,Dimy-1)):
            return True
        else:
            return False
        
    #This loop makes the random walker wiggle until the end condition
    while(CheckEnd(xPos,yPos)==False):
        xPos, yPos, xDistance, yDistance, SectorX,SectorY, lifetime = Advance(xPos,yPos,xDistance,yDistance,SectorX,SectorY,lifetime)
    return SectorX,SectorY, lifetime, xDistance, yDistance

#This loops over a handful of different lattice size.  The ** operate below is an exponential.
for j in xrange(6):   
    sectorx = 0
    sectory = 0
    secmix = 0
    totdist = 0
    lift = 0
    
    NumOfRuns = 1000000.
    
    for i in xrange(int(NumOfRuns)):
        
        secx,secy, t, xdisp,ydisp = SingleRun(8*2**j,8*2**j)
        if(secx==1 and secy==0):
            sectorx+=1
            totdist += math.sqrt(xdisp**2 + ydisp**2)
            lift+=t
        if(secy==1 and secx==0):
            sectory+=1
        if(secx==1 and secy == 1):
            secmix+=1
    
    print 8*2**j,sectorx / NumOfRuns,sectory / NumOfRuns, secmix / 100000., secmix / ((sectorx + sectory)/2.), totdist/(1.*sectorx),lift/(1.*sectorx)
