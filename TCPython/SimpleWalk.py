from scipy import random
import math

Directions = []

Directions.append((0,1))
Directions.append((0,-1))
Directions.append((1,0))
Directions.append((-1,0))





#Dimensions must be at least 5 by 5.
#Walker will start at (0,2).  If it becomes adjacent to (0,0), it will annihilate.  It will be assumed that the winding number
#axis is at least 2 columns away from the annihilation point.
def SingleRun(Dimx,Dimy):

    StartingPositions = []
    StartingPositions.append((0,2))
    StartingPositions.append((1,1))
    StartingPositions.append(((-1%Dimx),1))

    xPos,yPos = StartingPositions[random.randint(0,3)]
    xDistance,yDistance = 0,0
    lifetime = 0
    
    #xPos = 0
    #yPos = 2
    Sector = 0
    WindingAxis = math.floor(Dimx/2.)
    
    #print "what"
    
    def PrintState():
        for i in xrange(Dimx):
            print "\n"
            for j in xrange(Dimy):
                if (i,j) == (xPos,yPos):
                    print "1",
                else: print "0",
        print "\n"
    
    def Advance(xPos,yPos,xDistance,yDistance,Sector,lifetime):
        MoveX,MoveY = Directions[random.randint(0,4)]
        
        xDistance+=MoveX
        yDistance+=MoveY
        lifetime+=1
        
        if xPos == WindingAxis and MoveX == -1:
            Sector = (Sector + 1)%2
        if xPos == WindingAxis - 1 and MoveX == 1:
            Sector = (Sector + 1)%2
        
        xPos = (xPos + MoveX)%Dimx
        yPos = (yPos + MoveY)%Dimx
        
        return xPos,yPos,xDistance,yDistance,Sector,lifetime
    
    def CheckEnd(xPos,yPos):
        if ((xPos,yPos) == (1,0) or (xPos,yPos) == (0,1) or (xPos,yPos) == (Dimx-1,0) or (xPos,yPos) == (0,Dimy-1)):
            return True
        else:
            return False
        
    while(CheckEnd(xPos,yPos)==False):
        #print "blah"
        xPos, yPos, xDistance, yDistance, Sector, lifetime = Advance(xPos,yPos,xDistance,yDistance,Sector, lifetime)
        #PrintState()
    return Sector, lifetime, xDistance, yDistance

for j in [7,8]:   
    sector = 0
    totdist = 0
    lift = 0
    for i in xrange(100000):
        
        sec, t, xdisp,ydisp = SingleRun(8*2**j,8*2**j)
        if(sec==1):
            sector+=1
            totdist += math.sqrt(xdisp**2 + ydisp**2)
            lift+=t
    
    print 8*2**j,sector / 100000.,totdist/(1.*sector),lift/(1.*sector)
