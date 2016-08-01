import numpy
from copy import deepcopy
from numpy import *
from scipy import sparse
from scipy import random
#from numpy import linspace
import matplotlib.pyplot as plt
import time

import sys
import argparse as ap
import csv
from collections import deque

class Excitation:
    
    def __init__(self, xPos, yPos):
        self.x = xPos
        self.y = yPos
        
class Link:
    
    #convention will be that xPos1,yPos1 stores the vertex with the excitation (for translateable links)
    def __init__(self, xPos1, yPos1, xPos2, yPos2):
        self.x1 = xPos1
        self.y1 = yPos1
        self.x2 = xPos2
        self.y2 = yPos2
        
class Lattice:
    
    #You feed the lattice an array of 1s and 0s which will define locations of excitations
    #sector is the tuple [i,j] where i,j = {0,1}
    def __init__(self, ArrayOfOccupancies, sector):
        
        #Does not keep track of Createable links because randomly choosing a link on the lattice will probably result in an appropriate
        #choice for the regimes we're interested in.
        self.TranslateableList = []
        self.AnnihilateableList = []
        
        self.ExcitationList = []
        self.Sector = sector
        self.Dimension = len(ArrayOfOccupancies)
        
        for i in xrange(len(ArrayOfOccupancies)):
            for j in xrange(len(ArrayOfOccupancies)):
                if(ArrayOfOccupancies[i][j]==1):
                    NewExcitation = Excitation(i,j)
                    self.ExcitationList.append(NewExcitation)
                    
        for i in xrange(len(ArrayOfOccupancies)):
            for j in xrange(len(ArrayOfOccupancies)):
                '''if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[i][(j+1)%self.Dimension]==0): self.TranslateableList.append(Link(i,j,i,(j+1)%self.Dimension))
                if(ArrayOfOccupancies[i][j]==0 and ArrayOfOccupancies[i][(j+1)%self.Dimension]==1): self.TranslateableList.append(Link(i,(j+1)%self.Dimension,i,j))
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[(i+1)%self.Dimension][j]==0): self.TranslateableList.append(Link(i,j,(i+1)%self.Dimension,j))
                if(ArrayOfOccupancies[i][j]==0 and ArrayOfOccupancies[(i+1)%self.Dimension][j]==1): self.TranslateableList.append(Link((i+1)%self.Dimension,j,i,j))
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[i][(j+1)%self.Dimension]==1): self.AnnihilateableList.append(Link(i,j,i,(j+1)%self.Dimension))
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[(i+1)%self.Dimension][j]==1): self.AnnihilateableList.append(Link(i,j,(i+1)%self.Dimension))'''
                #vertex indices are NOT modded here so that the translation, creation, and annihilation methods can properly calculate the "distance" between adjacent
                #vertices (i.e. so they can perform sector checks)
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[i][(j+1)%self.Dimension]==0): self.TranslateableList.append(Link(i,j,i,(j+1)))
                if(ArrayOfOccupancies[i][j]==0 and ArrayOfOccupancies[i][(j+1)%self.Dimension]==1): self.TranslateableList.append(Link(i,(j+1),i,j))
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[(i+1)%self.Dimension][j]==0): self.TranslateableList.append(Link(i,j,(i+1),j))
                if(ArrayOfOccupancies[i][j]==0 and ArrayOfOccupancies[(i+1)%self.Dimension][j]==1): self.TranslateableList.append(Link((i+1),j,i,j))
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[i][(j+1)%self.Dimension]==1): self.AnnihilateableList.append(Link(i,j,i,(j+1)))
                if(ArrayOfOccupancies[i][j]==1 and ArrayOfOccupancies[(i+1)%self.Dimension][j]==1): self.AnnihilateableList.append(Link(i,j,(i+1)))
    
    def PrintState(self):
        myArray = numpy.zeros((self.Dimension,self.Dimension))
        
        for i in self.ExcitationList:
            myArray[i.x,i.y]+=1
            
        for i in xrange(len(myArray)):
            print "\n"
            for j in xrange(len(myArray)):
                print str(myArray[i,j]) + "\t",
        
        print "\n"
        print self.Sector
        print "\n"
        
    def ReturnOccupancies(self):
        myArray = numpy.zeros((self.Dimension,self.Dimension))
        
        for i in self.ExcitationList:
            myArray[i.x,i.y]+=1
        
        return myArray
        
        
    
    def SpawnExcitation(self, xPos, yPos):
        
        CollisionIndex = -1
        #If we're trying to add an excitation on top of one that already exists, annihilate the excitation
        for i in xrange(len(self.ExcitationList)):
            if self.ExcitationList[i].x == xPos and self.ExcitationList[i].y == yPos:
                CollisionIndex = i
                
        if(CollisionIndex!=-1):
            self.ExcitationList.pop(CollisionIndex)
        
        #Else, add the excitation
        if (CollisionIndex == -1):
            NewExcitation = Excitation(xPos,yPos)
            self.ExcitationList.append(NewExcitation)   
        
    def ReturnExcitationList(self):
        return self.ExcitationList
    
    def ReturnEnergy(self, Je):
        return 2*Je*len(self.ExcitationList)
    
    def ReturnExcitationIndex(self, XPos, YPos):
        for i in xrange(len(self.ExcitationList)):
            if self.ExcitationList[i].x == XPos and self.ExcitationList[i].y == YPos:
                return i
    
    #Given the coordinates of the excitation, checks to see if x1+dirx is greater or less than the dimension of the array
    #If it is, then a sector change occurs.  This function is multipurpose: if an excitation is being translated, then dirx and diry
    #represent the delta-x between initial and final locations.  For annihilation and creation, (dirx,diry) is just the offset to the other
    #excitation in the pair (to be created or destroyed).
    def CheckSector(self, x1,y1,dirx,diry):
        if x1+dirx == self.Dimension or \
        x1+dirx == -1:
            self.Sector[0] = (self.Sector[0] + 1) % 2
        
        if y1+diry == self.Dimension or \
        y1+diry == -1:
            self.Sector[1] = (self.Sector[1] + 1) % 2
    
    
    def TranslateExcitation(self, ExcitationIndex, XDisplacement, YDisplacement):
        ExcitationToMove = self.ExcitationList[ExcitationIndex]
        #Check to see if this movement will cause a change in Topological Sector (winding operators are
        #defined to look across the x=0 and y=0 axes (i.e., from x = Dimension to x = 0).
        self.CheckSector(ExcitationToMove.x,ExcitationToMove.y,XDisplacement,YDisplacement)
        
        ExcitationToMove.x = (ExcitationToMove.x + XDisplacement) % self.Dimension
        ExcitationToMove.y = (ExcitationToMove.y + YDisplacement) % self.Dimension
        
    
    def TranslateExcitations(self, ExcitationIndex, XDisplacement, YDisplacement):
        
        ExcitationToMove = self.ExcitationList[ExcitationIndex]
        #Check to see if this movement will cause a change in Topological Sector (winding operators are
        #defined to look across the x=0 and y=0 axes (i.e., from x = Dimension to x = 0).
        self.CheckSector(ExcitationToMove.x,ExcitationToMove.y,XDisplacement,YDisplacement)
        
        #print("**********")
        #print("Before:")
       # print(self.TranslateableList)
        
       # print("Deleting translateable links:")
        #print(self.TranslateableList)
        #print("**********")
        
        OldX = ExcitationToMove.x
        OldY = ExcitationToMove.y
        
        
        NewX = (ExcitationToMove.x + XDisplacement)%self.Dimension
        NewY = (ExcitationToMove.y + YDisplacement)%self.Dimension
            
        

	self.TranslateableList = [elem for elem in self.TranslateableList if not (elem.x1 == OldX and elem.y1 == OldY) and not (elem.x2 == NewX and elem.y2 == NewY)]

        self.AnnihilateableList = [elem for elem in self.AnnihilateableList if not ((elem.x1%self.Dimension == OldX and elem.y1%self.Dimension == OldY) \
                                                                                    or (elem.x2%self.Dimension == OldX and elem.y2%self.Dimension == OldY))]
        
        
        ExcitationToMove.x = (ExcitationToMove.x + XDisplacement) % self.Dimension
        ExcitationToMove.y = (ExcitationToMove.y + YDisplacement) % self.Dimension
        
        Occupancies = self.ReturnOccupancies()
        
        if Occupancies[(OldX+1)%self.Dimension,OldY] == 1 and ((OldX+1)%self.Dimension,OldY) != (NewX,NewY): self.TranslateableList.append(Link((OldX+1)%self.Dimension,OldY,(OldX+1)%self.Dimension-1,OldY))
        if Occupancies[OldX,(OldY+1)%self.Dimension] == 1 and (OldX,(OldY+1)%self.Dimension) != (NewX,NewY): self.TranslateableList.append(Link(OldX,(OldY+1)%self.Dimension,OldX,(OldY+1)%self.Dimension-1))
        if Occupancies[(OldX-1)%self.Dimension,OldY] == 1 and ((OldX-1)%self.Dimension,OldY) != (NewX,NewY): self.TranslateableList.append(Link((OldX-1)%self.Dimension,OldY,(OldX-1)%self.Dimension+1,OldY))
        if Occupancies[OldX,(OldY-1)%self.Dimension] == 1 and (OldX,(OldY-1)%self.Dimension) != (NewX,NewY): self.TranslateableList.append(Link(OldX,(OldY-1)%self.Dimension,OldX,(OldY-1)%self.Dimension+1))
        
        if Occupancies[NewX,(NewY+1)%self.Dimension]==0:
            self.TranslateableList.append(Link(NewX,NewY,NewX,NewY+1))
        else:
            self.AnnihilateableList.append(Link(NewX,NewY,NewX,NewY+1))
            
        if Occupancies[(NewX+1)%self.Dimension,NewY]==0:
            self.TranslateableList.append(Link(NewX,NewY,NewX+1,NewY))
        else:
            self.AnnihilateableList.append(Link(NewX,NewY,NewX+1,NewY))
            
        if Occupancies[NewX,(NewY-1)%self.Dimension]==0:
            self.TranslateableList.append(Link(NewX,NewY,NewX,NewY-1))
        else:
            self.AnnihilateableList.append(Link(NewX,NewY,NewX,NewY-1))
            
        if Occupancies[(NewX-1)%self.Dimension,NewY]==0:
            self.TranslateableList.append(Link(NewX,NewY,NewX-1,NewY))
        else:
            self.AnnihilateableList.append(Link(NewX,NewY,NewX-1,NewY))
        
        
        #print("**********")
        #print("Populated translateable links:")
        #print(self.TranslateableList)
        #print("**********")
        
        
    
    
    #These coordinate are NOT modulo the dimension.  This is so that the sectorcheck can work properly.
    def AnnihilateExcitations(self,xPos1,yPos1,xPos2,yPos2):
        
        self.AnnihilateableList = [elem for elem in self.AnnihilateableList if not (((elem.x1%self.Dimension==xPos1%self.Dimension and elem.y1%self.Dimension==yPos1%self.Dimension) \
                                   or (elem.x1%self.Dimension==xPos2%self.Dimension and elem.y1%self.Dimension==yPos2%self.Dimension)) \
                                or ((elem.x2%self.Dimension==xPos1%self.Dimension and elem.y2%self.Dimension==yPos1%self.Dimension) \
                                   or (elem.x2%self.Dimension==xPos2%self.Dimension and elem.y2%self.Dimension==yPos2%self.Dimension)))]
        self.TranslateableList = [elem for elem in self.TranslateableList if not ((elem.x1 == xPos1%self.Dimension and elem.y1 == yPos1%self.Dimension)\
                                  or (elem.x1 == xPos2%self.Dimension and elem.y1 == yPos2%self.Dimension))]
        
        
        '''print("**********")
        print("After annihilation, translateable links:")
        print(self.TranslateableList)
        print("**********")
    
        print("**********")
        print("After annihilation, annihilateable links:")
        print(self.AnnihilateableList)
        print("**********")'''
    
    
        self.CheckSector(xPos1%self.Dimension,yPos1%self.Dimension,xPos2-xPos1,yPos2-yPos1)
        pop1 = self.ReturnExcitationIndex(xPos1%self.Dimension,yPos1%self.Dimension)
        pop2 = self.ReturnExcitationIndex(xPos2%self.Dimension,yPos2%self.Dimension)
        
        self.ExcitationList.pop(max(pop1,pop2))
        self.ExcitationList.pop(min(pop1,pop2))
        
    def CreateExcitations(self,xPos1,yPos1,xPos2,yPos2):
        
        
        
        self.CheckSector(xPos1%self.Dimension,yPos1%self.Dimension,xPos2-xPos1,yPos2-yPos1)
        self.ExcitationList.append(Excitation(xPos1%self.Dimension,yPos1%self.Dimension))
        self.ExcitationList.append(Excitation(xPos2%self.Dimension,yPos2%self.Dimension))
        
        self.GenerateAnnihilationLinks(xPos1%self.Dimension,yPos1%self.Dimension)
        self.GenerateAnnihilationLinks(xPos2%self.Dimension,yPos2%self.Dimension)
        #account for doublecounting
        self.AnnihilateableList = [elem for elem in self.AnnihilateableList if not ((elem.x1%self.Dimension==xPos1%self.Dimension \
                                                                                and elem.y1%self.Dimension==yPos1%self.Dimension\
                                                                               and elem.x2%self.Dimension==xPos2%self.Dimension \
                                                                               and elem.y2%self.Dimension==yPos2%self.Dimension) or\
                                                                               (elem.x2%self.Dimension==xPos1%self.Dimension \
                                                                               and elem.y2%self.Dimension==yPos1%self.Dimension\
                                                                               and elem.x1%self.Dimension==xPos2%self.Dimension \
                                                                               and elem.y1%self.Dimension==yPos2%self.Dimension))]
        self.AnnihilateableList.append(Link(xPos1,yPos1,xPos2,yPos2))
        
                              
        
        self.GenerateTranslationLinks(xPos1%self.Dimension,yPos1%self.Dimension)
        self.GenerateTranslationLinks(xPos2%self.Dimension,yPos2%self.Dimension)
        
    def AnnihilateExcitation(self, Xpos, Ypos):
        self.ExcitationList.pop(self.ReturnExcitationIndex(Xpos,Ypos))
        
    def GenerateAnnihilationLinks(self,xPos,yPos):
        Occupancies = self.ReturnOccupancies()
        
        if Occupancies[xPos,(yPos+1)%self.Dimension]==1:
            self.AnnihilateableList.append(Link(xPos,yPos,xPos,yPos+1))
            
        if Occupancies[(xPos+1)%self.Dimension,yPos]==1:
            self.AnnihilateableList.append(Link(xPos,yPos,xPos+1,yPos))
            
        if Occupancies[xPos,(yPos-1)%self.Dimension]==1:
            self.AnnihilateableList.append(Link(xPos,yPos,xPos,yPos-1))
            
        if Occupancies[(xPos-1)%self.Dimension,yPos]==1:
            self.AnnihilateableList.append(Link(xPos,yPos,xPos-1,yPos))
    
    def GenerateTranslationLinks(self,xPos,yPos):
        Occupancies = self.ReturnOccupancies()
        
        if Occupancies[xPos,(yPos+1)%self.Dimension]==0:
            self.TranslateableList.append(Link(xPos,yPos,xPos,yPos+1))
            
        if Occupancies[(xPos+1)%self.Dimension,yPos]==0:
            self.TranslateableList.append(Link(xPos,yPos,xPos+1,yPos))
            
        if Occupancies[xPos,(yPos-1)%self.Dimension]==0:
            self.TranslateableList.append(Link(xPos,yPos,xPos,yPos-1))
            
        if Occupancies[(xPos-1)%self.Dimension,yPos]==0:
            self.TranslateableList.append(Link(xPos,yPos,xPos-1,yPos))
        
        
        

    
    
    
    
    #indexinlist is just the local index of the excitation we want to move, as stored in ExcitationList
    #direction is 0,1,2,3, and is just North, East, South, or West with (0,0) in the bottom left corner
    #returns True if no collision occurred, else false
    def MoveExcitation(self, IndexInList, direction):
        
        #OrderedPair stores the actual change in position in (x,y) space (i.e., a mapping of direction to [+x,+y]
        OrderedPair = [0,0]
        if direction == 0:
            OrderedPair = [0,1]
        if direction == 1:
            OrderedPair = [1,0]
        if direction == 2:
            OrderedPair = [0,-1]
        if direction == 3:
            OrderedPair = [-1,0]
        
        ExcitationToMove = self.ExcitationList[IndexInList]    
        
        #Check to see if this movement will cause a change in Topological Sector (winding operators are
        #defined to look across the x=0 and y=0 axes (i.e., from x = Dimension to x = 0).
        self.CheckSector(ExcitationToMove.x,ExcitationToMove.y,ExcitationToMove.x+OrderedPair[0],ExcitationToMove.y+OrderedPair[1])
        
        
        
        #Check to see if another excitation is already where we want to move
        Collision = -1
        
        for i in xrange(len(self.ExcitationList)):
            #Checks to see if we're going to move our excitation on top of another excitation
            if (i!=IndexInList and self.ExcitationList[i].x == (ExcitationToMove.x + OrderedPair[0])%self.Dimension \
            and self.ExcitationList[i].y == (ExcitationToMove.y + OrderedPair[1])%self.Dimension):
                Collision = i
        
        #If there is an excitation there, annihilate the excitations
        if Collision != -1:
            self.ExcitationList.pop(max(Collision, IndexInList))
            self.ExcitationList.pop(min(Collision, IndexInList))
            return False
        else:
            ExcitationToMove.x = (ExcitationToMove.x + OrderedPair[0]) % self.Dimension
            ExcitationToMove.y = (ExcitationToMove.y + OrderedPair[1]) % self.Dimension
            return True
            

def SimulateSingleShot(MyLattice, SpawnRate, AnnihilationRate, TranslateRate, CurrentTime, TotalTime):
    Times = []
    WindingObservables = []
    Occupations = [2]
    Times.append(CurrentTime)
    WindingObservables.append(list(MyLattice.Sector))
    
    def AdvanceTime(TimeList,ObservableList):
        
        TranslateableLinks = MyLattice.TranslateableList
        AnnihilateableLinks = MyLattice.AnnihilateableList
        
        Norm = len(TranslateableLinks)*TranslateRate + len(AnnihilateableLinks)*AnnihilationRate + \
        ((2*(MyLattice.Dimension)**2) - len(AnnihilateableLinks) - len(TranslateableLinks))*SpawnRate
        
        PTrans = len(TranslateableLinks)*TranslateRate/Norm
        PAnnih = len(AnnihilateableLinks)*AnnihilationRate/Norm
        PSpawn = ((2*(MyLattice.Dimension)**2) - len(AnnihilateableLinks) - len(TranslateableLinks))*SpawnRate/Norm
        
        
        r = random.rand()
        DeltaTau = (-1./Norm)*log(r)
        
        if (r < PTrans):
            RandLink = random.randint(len(TranslateableLinks))
            LinkCoords = [TranslateableLinks[RandLink].x1,TranslateableLinks[RandLink].y1,TranslateableLinks[RandLink].x2,TranslateableLinks[RandLink].y2]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            MyLattice.TranslateExcitations(MyLattice.ReturnExcitationIndex(LinkCoords[0],LinkCoords[1]), XDisplacement, YDisplacement)
            
        if (r >= PTrans and r < PTrans + PAnnih):
            RandLink = random.randint(len(AnnihilateableLinks))
            LinkCoords = [AnnihilateableLinks[RandLink].x1,AnnihilateableLinks[RandLink].y1,AnnihilateableLinks[RandLink].x2,AnnihilateableLinks[RandLink].y2]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            MyLattice.AnnihilateExcitations(LinkCoords[0],LinkCoords[1],LinkCoords[2],LinkCoords[3])

        if (r >= PTrans + PAnnih):
            Occupancies = MyLattice.ReturnOccupancies()
            FoundCreatableLink = False
            RandX = 0
            RandY = 0
            RandDir = 0
            
            while FoundCreatableLink == False:
                RandX = random.randint(MyLattice.Dimension)
                RandY = random.randint(MyLattice.Dimension)
                RandDir = random.randint(2)
                if(Occupancies[RandX,RandY]==0 and \
                   Occupancies[(RandX + RandDir)%MyLattice.Dimension,(RandY + (1-RandDir))%MyLattice.Dimension] == 0):
                    FoundCreatableLink = True
        
            LinkCoords = [RandX,RandY,RandX + RandDir, RandY + (1-RandDir)]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]

            MyLattice.CreateExcitations(LinkCoords[0],LinkCoords[1],LinkCoords[2],LinkCoords[3])
        
        return DeltaTau, MyLattice.Sector, len(MyLattice.ExcitationList)
    
    while Occupations[-1] == 2:
        TimeIncrement, Sector, NumOfExcitations = AdvanceTime(Times,WindingObservables)
        Times.append(Times[-1] + TimeIncrement)
        WindingObservables.append(list(Sector))
        Occupations.append(NumOfExcitations)
        #MyLattice.PrintState()
        #print(Times[-1])

    return Times, WindingObservables, Occupations




def Simulate3(MyLattice, SpawnRate, AnnihilationRate, TranslateRate, CurrentTime, TotalTime):
    Times = []
    WindingObservables = []
    Occupations = [0]
    Times.append(CurrentTime)
    WindingObservables.append(list(MyLattice.Sector))
    
    
    def AdvanceTime(TimeList,ObservableList):
        
        #Norm = len(Trans)*TranslateRate + len(Annih)*AnnihilationRate + ((2*(MyLattice.Dimension)**2) - len(Annih) - len(Trans))*SpawnRate
        
        
        
        #print "Advancing time!"
        
        TranslateableLinks = MyLattice.TranslateableList
        AnnihilateableLinks = MyLattice.AnnihilateableList
        
        #print "Translateable: " + str(len(TranslateableLinks)) + str(TranslateableLinks)
        #print "Annihilateable: " + str(len(AnnihilateableLinks)) + str(AnnihilateableLinks)
        
        Norm = len(TranslateableLinks)*TranslateRate + len(AnnihilateableLinks)*AnnihilationRate + \
        ((2*(MyLattice.Dimension)**2) - len(AnnihilateableLinks) - len(TranslateableLinks))*SpawnRate
        
        PTrans = len(TranslateableLinks)*TranslateRate/Norm
        PAnnih = len(AnnihilateableLinks)*AnnihilationRate/Norm
        PSpawn = ((2*(MyLattice.Dimension)**2) - len(AnnihilateableLinks) - len(TranslateableLinks))*SpawnRate/Norm
        
        
        r = random.rand()
        DeltaTau = (-1./Norm)*log(r)
        
        if (r < PTrans):
            
            #print "Translating!"
            #MyLattice.PrintState()
            RandLink = random.randint(len(TranslateableLinks))
            LinkCoords = [TranslateableLinks[RandLink].x1,TranslateableLinks[RandLink].y1,TranslateableLinks[RandLink].x2,TranslateableLinks[RandLink].y2]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            MyLattice.TranslateExcitations(MyLattice.ReturnExcitationIndex(LinkCoords[0],LinkCoords[1]), XDisplacement, YDisplacement)
            
        if (r >= PTrans and r < PTrans + PAnnih):
            
            #print "Annihilating!"
            
            RandLink = random.randint(len(AnnihilateableLinks))
            LinkCoords = [AnnihilateableLinks[RandLink].x1,AnnihilateableLinks[RandLink].y1,AnnihilateableLinks[RandLink].x2,AnnihilateableLinks[RandLink].y2]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            #if(LinkCoords[2]==MyLattice.Dimension or LinkCoords[0] == MyLattice.Dimension - 1):
                #print str(LinkCoords[2]) + "\t" + str(LinkCoords[0])
                #print str(LinkCoords[3]) + "\t" + str(LinkCoords[1])
            MyLattice.AnnihilateExcitations(LinkCoords[0],LinkCoords[1],LinkCoords[2],LinkCoords[3])
            #MyLattice.CheckSector(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension, XDisplacement, YDisplacement)
        
        if (r >= PTrans + PAnnih):
            
            #print "Creating!"
            
            Occupancies = MyLattice.ReturnOccupancies()
            FoundCreatableLink = False
            RandX = 0
            RandY = 0
            RandDir = 0
            
            while FoundCreatableLink == False:
                RandX = random.randint(MyLattice.Dimension)
                RandY = random.randint(MyLattice.Dimension)
                RandDir = random.randint(2)
                if(Occupancies[RandX,RandY]==0 and \
                   Occupancies[(RandX + RandDir)%MyLattice.Dimension,(RandY + (1-RandDir))%MyLattice.Dimension] == 0):
                    FoundCreatableLink = True
        
            LinkCoords = [RandX,RandY,RandX + RandDir, RandY + (1-RandDir)]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            #if(LinkCoords[2]==MyLattice.Dimension or LinkCoords[0] == MyLattice.Dimension - 1):
                #print str(LinkCoords[2]) + "\t" + str(LinkCoords[0])
                #print str(LinkCoords[3]) + "\t" + str(LinkCoords[1])
            MyLattice.CreateExcitations(LinkCoords[0],LinkCoords[1],LinkCoords[2],LinkCoords[3])
            #MyLattice.CheckSector(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension, XDisplacement, YDisplacement)
        
        return DeltaTau, MyLattice.Sector, len(MyLattice.ExcitationList)
    
    '''OutputMode = "Slim"
    
    if OutputMode == "Verbose":
        while Times[-1] < TotalTime:
            TimeIncrement, Sector, NumOfExcitations = AdvanceTime(Times,WindingObservables)
            Times.append(Times[-1] + TimeIncrement)
            WindingObservables.append(list(Sector))
            Occupations.append(NumOfExcitations)
            #print([Times[-1],WindingObservables[-1]])
            #MyLattice.PrintState()
            #print(Times[-1])
        #print Occupations
            
    if OutputMode == "Slim":   
        TimePartition = numpy.linspace(0,TotalTime,1000)
        #print TimePartition
        TimeCounter = 0
        Delay = 1
        while CurrentTime < TotalTime:
            TimeIncrement, Sector, NumOfExcitations = AdvanceTime(Times,WindingObservables)
            CurrentTime = CurrentTime + TimeIncrement
            
            if CurrentTime >= TimePartition[TimeCounter] and TimeCounter < 998:
                if Delay == 0:
                #Times.append(Times[-1] + TimeIncrement)
                    Times.append(float(CurrentTime))
                    WindingObservables.append(list(Sector))
                    Occupations.append(NumOfExcitations)
                    #print NumOfExcitations
                    TimeCounter+=1
                    Delay+=2
                Delay = Delay - 1'''    

    while Times[-1] < TotalTime:
        TimeIncrement, Sector, NumOfExcitations = AdvanceTime(Times,WindingObservables)
        Times.append(Times[-1] + TimeIncrement)
        WindingObservables.append(list(Sector))
        Occupations.append(NumOfExcitations)
        #MyLattice.PrintState()
        #print(Times[-1])

    return Times, WindingObservables, Occupations




def Simulate2(MyLattice, SpawnRate, AnnihilationRate, TranslateRate, CurrentTime, TotalTime):
    
    Times = []
    Observables = []
    Times.append(CurrentTime)
    Observables.append(list(MyLattice.Sector))
    
    
    def AdvanceTime(TimeList,ObservableList):
        
        PTrans, PAnnih, PSpawn, Norm, TranslateableLinks, AnnihilateableLinks, SpawnableLinks = Probabilities()
        
        r = random.rand()
        DeltaTau = (-1./Norm)*log(r)
        
        if (r < PTrans):
            RandLink = random.randint(len(TranslateableLinks))
            LinkCoords = TranslateableLinks[RandLink]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            MyLattice.TranslateExcitation(MyLattice.ReturnExcitationIndex(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension), XDisplacement, YDisplacement)
            
        if (r >= PTrans and r < PTrans + PAnnih):
            RandLink = random.randint(len(AnnihilateableLinks))
            LinkCoords = AnnihilateableLinks[RandLink]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            #if(LinkCoords[2]==MyLattice.Dimension or LinkCoords[0] == MyLattice.Dimension - 1):
                #print str(LinkCoords[2]) + "\t" + str(LinkCoords[0])
                #print str(LinkCoords[3]) + "\t" + str(LinkCoords[1])
            MyLattice.AnnihilateExcitation(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension)
            MyLattice.AnnihilateExcitation(LinkCoords[2]%MyLattice.Dimension,LinkCoords[3]%MyLattice.Dimension)
            MyLattice.CheckSector(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension, XDisplacement, YDisplacement)
        
        if (r >= PTrans + PAnnih):
            RandLink = random.randint(len(SpawnableLinks))
            LinkCoords = SpawnableLinks[RandLink]
            XDisplacement = LinkCoords[2] - LinkCoords[0]
            YDisplacement = LinkCoords[3] - LinkCoords[1]
            #if(LinkCoords[2]==MyLattice.Dimension or LinkCoords[0] == MyLattice.Dimension - 1):
                #print str(LinkCoords[2]) + "\t" + str(LinkCoords[0])
                #print str(LinkCoords[3]) + "\t" + str(LinkCoords[1])
            MyLattice.SpawnExcitation(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension)
            MyLattice.SpawnExcitation(LinkCoords[2]%MyLattice.Dimension,LinkCoords[3]%MyLattice.Dimension)
            MyLattice.CheckSector(LinkCoords[0]%MyLattice.Dimension,LinkCoords[1]%MyLattice.Dimension, XDisplacement, YDisplacement)
        
        return DeltaTau, MyLattice.Sector
        
        
        #Returns tuples for the links that could be translated or annihilated
        #Finding links for creation is just a matter of picking a pair that isn't contained
        #in this list.
    
        #Convention is that [x,y,x',y'] is defined such that x,y is the vertex with an excitation
        #(for translateable links, obviously)
    def LinkLocations(TypeOfLink):
            
        LinksToReturn = []
        
        Occupancies = MyLattice.ReturnOccupancies()
        
        for x in xrange(MyLattice.Dimension):
            for y in xrange(MyLattice.Dimension):
                
                #print Translatables
                
                CommonVertex = Occupancies[x,y]
                AboveVertex = Occupancies[x,(y+1)%MyLattice.Dimension]
                RightVertex = Occupancies[(x+1)%MyLattice.Dimension,y]
                
                #The translatable coordinates aren't modded because what's important are differences
                #between coordinates.  If the modding was done here, bad things would happen, for example:
                #in a, say, 8x8 lattice, [0,7] and [0,0] would be 7 apart.
                #The modding is done above when the excitation is actually being added
                if(TypeOfLink == "Translator"):
                    if(CommonVertex != AboveVertex):
                        if(CommonVertex == 1):
                            LinksToReturn.append([x,y,x,(y+1)])
                        else:
                            LinksToReturn.append([x,(y+1),x,y])
                        
                    if(CommonVertex != RightVertex):
                        if(CommonVertex == 1):
                            LinksToReturn.append([x,y,(x+1),y])
                        else:
                            LinksToReturn.append([(x+1),y,x,y])
                    
                if(TypeOfLink == "Annihilator"):
                    
                    if(CommonVertex == 1 and AboveVertex == CommonVertex):
                        LinksToReturn.append([x,y,x,(y+1)])
                        
                    if(CommonVertex == 1 and RightVertex == CommonVertex):
                        LinksToReturn.append([x,y,(x+1),y])
                
                if(TypeOfLink == "Creator"):
                    if(CommonVertex == 0 and AboveVertex == 0):
                        LinksToReturn.append([x,y,x,(y+1)])
                    if(CommonVertex == 0 and RightVertex == 0):
                        LinksToReturn.append([x,y,(x+1),y])
        
        return LinksToReturn
        
        
    #Returns the tuple (prob of translation, prob of annihilation, prob of spawning, Norm)
    def Probabilities():
        
        Trans = LinkLocations("Translator")
        Annih = LinkLocations("Annihilator")
        Spawnables = LinkLocations("Creator")
        Norm = len(Trans)*TranslateRate + len(Annih)*AnnihilationRate + ((2*(MyLattice.Dimension)**2) - len(Annih) - len(Trans))*SpawnRate
        
        return len(Trans)*TranslateRate/Norm, len(Annih)*AnnihilationRate/Norm, ((2*(MyLattice.Dimension)**2) - len(Annih) - len(Trans))*SpawnRate/Norm, Norm, Trans, Annih, Spawnables
    

    while Times[-1] < TotalTime:
        TimeIncrement, Sector = AdvanceTime(Times,Observables)
        Times.append(Times[-1] + TimeIncrement)
        Observables.append(list(Sector))
        #MyLattice.PrintState()
        #print(Times[-1])

    return Times, Observables, _



#Je is the term in the toric code hamiltonian for the vertex operator (governing electric excitations
def Simulate(MyLattice, Je, Temperature, Timesteps):
    
    #At each time step, each pair of adjacent coordinates could create an excitation, and
    #each existing excitation could move (or not move).  We're essentially implementing the 
    #metropolis algorithm for nearby configurations.

    def TrialMove2(InitLat):
        InitialLattice = deepcopy(InitLat)
        InitialEnergy = 2*Je*len(InitialLattice.ExcitationList)
        
        if(len(InitialLattice.ExcitationList)!=0):
            RandomExcitation = random.randint(len(InitialLattice.ExcitationList))
            RandomDirection = random.randint(0,5)
            
            InitialLattice.MoveExcitation(RandomExcitation,RandomDirection)
            
        SpawnOrientation = random.randint(0,3) #0 don't do anything, 1 updown, 2 leftright
        
        if(SpawnOrientation == 1):
            RandomX = random.randint(0,InitialLattice.Dimension)
            RandomY = random.randint(0,InitialLattice.Dimension)
            
            InitialLattice.SpawnExcitation(RandomX,RandomY)
            InitialLattice.SpawnExcitation(RandomX,(RandomY+1)%InitialLattice.Dimension)
            InitialLattice.CheckSector(RandomX,RandomY,0,1)
            
        if(SpawnOrientation == 2):
            RandomX = random.randint(0,InitialLattice.Dimension)
            RandomY = random.randint(0,InitialLattice.Dimension)
            
            InitialLattice.SpawnExcitation(RandomX,RandomY)
            InitialLattice.SpawnExcitation((RandomX+1)%InitialLattice.Dimension,RandomY)
            InitialLattice.CheckSector(RandomX,RandomY,1,0)
            
        FinalEnergy = 2*Je*len(InitialLattice.ExcitationList)
        
        if(FinalEnergy <= InitialEnergy):
            return InitialLattice
        
        else:
            r = random.rand()
            #print(exp(-(FinalEnergy-InitialEnergy)/Temperature))
            if r < exp(-(FinalEnergy-InitialEnergy)/Temperature):
                return InitialLattice
            else:
                return InitLat

        
        
    def TrialMove(InitLat):
        InitialLattice = deepcopy(InitLat)
        InitialEnergy = 2*Je*len(InitialLattice.ExcitationList)
        
        i = 0
        ListLength = len(InitialLattice.ExcitationList)
        while i < ListLength:
            RandomDirection = random.randint(0,5)
            if InitialLattice.MoveExcitation(i,RandomDirection) == False:
                #if there was a collision, and this collision killed all the excitations
                #exit the loop
                if len(InitialLattice.ExcitationList) == 0:
                    i = 1
                #else, you've just popped the excitation under consideration, so now i is
                #point to something else
                else:
                    i = i
            #if you didn't collide with anything, continue iterating over excitations
            else:
                i += 1
            
            #The number of excitations may have gone down, so we need to keep track of that:
            ListLength = len(InitialLattice.ExcitationList)
        
        for x in xrange(InitialLattice.Dimension):
            for y in xrange(InitialLattice.Dimension):
                SpawnPairUpDown = random.randint(0,2)
                SpawnPairLeftRight = random.randint(0,2)
                
                if SpawnPairUpDown == 1:
                    InitialLattice.SpawnExcitation(x,y)
                    InitialLattice.SpawnExcitation(x,(y+1)%InitialLattice.Dimension)
                    InitialLattice.CheckSector(x,y,0,1)
                
                if SpawnPairLeftRight == 1:
                    InitialLattice.SpawnExcitation(x,y)
                    InitialLattice.SpawnExcitation((x+1)%InitialLattice.Dimension,y)
                    InitialLattice.CheckSector(x,y,1,0)
                    
        FinalEnergy = 2*Je*len(InitialLattice.ExcitationList)
        
        if(FinalEnergy <= InitialEnergy):
            return InitialLattice
        
        else:
            r = random.rand()
            print(exp(-(FinalEnergy-InitialEnergy)/Temperature))
            if r < exp(-(FinalEnergy-InitialEnergy)/Temperature):
                return InitialLattice
            else:
                return InitLat

        
    TimeInSector00 = []
    TimeInSector01 = []
    TimeInSector10 = []
    TimeInSector11 = []     
    
    counter = 0
    
    for t in xrange(Timesteps):
        if t%10000 == 0:
            print t
        SectorBefore = MyLattice.Sector
        MyLattice = TrialMove2(MyLattice)
        SectorAfter = MyLattice.Sector
        
        if(SectorBefore == SectorAfter):
            counter+=1
        else:
            if SectorBefore == [0,0]: TimeInSector00.append(counter)
            if SectorBefore == [0,1]: TimeInSector01.append(counter)
            if SectorBefore == [1,0]: TimeInSector10.append(counter)
            if SectorBefore == [1,1]: TimeInSector11.append(counter)
            counter = 0
        
        #print(MyLattice.ReturnEnergy(Je))
        #MyLattice.PrintState()
        
        if(t==Timesteps - 1):
            if SectorBefore == [0,0]: TimeInSector00.append(counter)
            if SectorBefore == [0,1]: TimeInSector01.append(counter)
            if SectorBefore == [1,0]: TimeInSector10.append(counter)
            if SectorBefore == [1,1]: TimeInSector11.append(counter)
        
    print(TimeInSector00)
    print("***********")
    print(TimeInSector01)
    print("***********")
    print(TimeInSector10)
    print("***********")
    print(TimeInSector11)
    print("***********")
    #print(sum(TimeInSector00)/len(TimeInSector00)) 
    #print(sum(TimeInSector01)/len(TimeInSector01))
    #print(sum(TimeInSector10)/len(TimeInSector10))
    #print(sum(TimeInSector11)/len(TimeInSector11))
    #print(str(sum(TimeInSector00) + sum(TimeInSector10) + sum(TimeInSector01) + sum(TimeInSector11)))
    
    return [TimeInSector00,TimeInSector01,TimeInSector10,TimeInSector11]


    
    #print CurrentIndexes
    #print SetOfTimesAndObservables[0][0]
    #print SetOfTimesAndObservables[1][0]

#for t in xrange(len(DataList)):
#    print str(DataList[t][0]) + "\t" + str(DataList[t][1])
                                 

    



'''for t in Temperatures:
    ThisTempAverage = []
    ThisTempStDev = []
    MyLattice = Lattice(ListOfExcitations,[0,0])
    TimesInSectors = Simulate(MyLattice, 1, t, 10000000) #10000000
    for i in xrange(len(TimesInSectors)):
        if len(TimesInSectors[i])!=0:
            ThisAverage = sum(TimesInSectors[i])/len(TimesInSectors[i])
            ThisTempAverage.append(ThisAverage)
            StdDev = 0
            for j in xrange(len(TimesInSectors[i])):
                StdDev += (TimesInSectors[i][j])**2
            StdDev = sqrt((StdDev / len(TimesInSectors[i])) - ThisAverage**2)
            ThisTempStDev.append(StdDev)
        else:
            ThisTempAverage.append(0)
            ThisTempStDev.append(0)
    Averages.append(ThisTempAverage)
    Deviations.append(ThisTempStDev)
    
plt.clf()
plt.plot(Temperatures, Averages)
plt.xlabel(r"Temperature")
plt.ylabel(r"Time in Sector")
plt.ylim((0,2000000))
plt.title(r"Expectation Values in each Sector at given Temperature")
plt.savefig("plots/averages" + str(time.time()) + ".png", format="png")

plt.clf()
plt.plot(Temperatures, Deviations)
plt.xlabel(r"Temperature")
plt.ylabel(r"Standard Deviation of Time")
plt.ylim((0,2000000))
plt.title(r"Standard Deviations of Expected Time In Each Sector")
plt.savefig("plots/stdevs" + str(time.time()) + ".png", format="png")

for i in xrange(len(Averages)):
    print "\n"
    for j in xrange(len(Averages[i])):
        print(str(Averages[i][j]) + "\t"),

print "\n"
print("**********************")

print Temperatures
print Averages
print Deviations

with open("averages" + ".csv", "w") as outFile:
    writer = csv.writer(outFile)
    for Temp, Avg in zip(Temperatures, Averages):
        row = deque(Avg)
        row.appendleft(Temp)
        writer.writerow(row)
        
with open("deviations" + ".csv", "w") as outFile:
    writer = csv.writer(outFile)
    for Temp, Dev in zip(Temperatures, Deviations):
        row = deque(Dev)
        row.appendleft(Temp)
        writer.writerow(row)'''
    
            
    
    
