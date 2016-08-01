"""
Ising model one-shot dynamics simulation.
From C. Daniel Freeman (2016 http://arxiv.org/abs/1603.05005)
"""

import logging
import math
import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
#import isingutils.py

import random
from random import choice
import copy
import sys
from compiler.ast import flatten
from numpy import *



logger = logging.getLogger(__name__)

class IsingEnv(gym.Env):
	metadata = {
		'render.modes': ['human', 'rgb_array'],
		'video.frames_per_second' : 50
	}
	
	def __init__(self):
	
		#Holds transform objects for rendering
		self.translist = []
		self.error_translist = []
	
		self.TotalTime = 0.
		#self.NextActionTime = 0.
		#self.NextActionProbability = 0.
		
		self.SystemLength = 24
		self.Temperature = .15
		self.Delta = 1.0
		self.CreationRate = abs(1./(1-np.exp(self.Delta*1.0/self.Temperature)))
		self.AnnihilationRate = abs(1./(1-np.exp(-self.Delta*1.0/self.Temperature)))
		self.HoppingRate = .01#self.Temperature
		self.CorrectionRate = 1.
		
		self.Sector = 0
		self.state = np.zeros(self.SystemLength)
		self.error_state = np.zeros(self.SystemLength)

		
		# Angle limit set to 2 * theta_threshold_radians so failing observation is still within bounds
		low = np.zeros(self.SystemLength)
		high = np.ones(self.SystemLength)
		# Can perform a swap between any pair of sites.  Convention is that 0 swaps from 0 to 1 and (SystemLength-1) swaps from (SystemLength-1) to 0.
		# i.e., periodic boundary conditions.
		self.action_space = spaces.Discrete(self.SystemLength)
		self.observation_space = spaces.Box(low, high)
		self._seed()
		self.reset()
		self.viewer = None
		
		anyons_list = self.state

		self.steps_beyond_done = 0. 
		
		#Need to calculate when the first bath interaction will occur
		ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(anyons_list)
		Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate + len(ExPairLocList)*self.AnnihilationRate + \
		(len(EmptyPairLocList))*self.CreationRate
		
		self.PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate/Norm
		self.PAnn = len(ExPairLocList)*self.AnnihilationRate/Norm
		self.PCre = (len(anyons_list) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*self.HoppingRate/Norm
				
		self.NextActionProbability = random.random()
		self.NextActionTime = self.TotalTime + (-1./Norm)*np.log(self.NextActionProbability)

		# Just need to initialize the relevant attributes
		self._configure()

	def _configure(self, display=None):
		self.display = display

	def _seed(self, seed=None):
		self.np_random, seed = seeding.np_random(seed)
		return [seed]

	def _step(self, action):
		assert self.action_space.contains(action), "%r (%s) invalid"%(action, type(action))
		
		
		print "Current time: " + str(self.TotalTime)
		print "Next action: " + str(self.NextActionTime)
		state = self.state
		
		anyons_list = state
		
		#Store the winding operator before we do anything to the chain
		p = int(floor(len(self.state)/2.))
		pl,pr = self.state[p],self.state[p+1] 
		
		#I'm going to change this into a more discrete picture--where there's an integer system clock,
		#and dynamics occur inbetween clock calls.
		
		#The most obvious way to do this is to stick a while loop in the if statement below that performs dynamics until
		#the next action time is after the next cycle of (CorrectionPeriod) * (n + 1)  (if this were occuring between n and n+1)
		
		#if the next corrective action would occur after the next bath interaction, do the bath interaction and calculate the next bath interaction time
		if self.TotalTime + (1. / self.CorrectionRate) > self.NextActionTime:
			ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(self.state)
			Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate + len(ExPairLocList)*self.AnnihilationRate + \
			(len(EmptyPairLocList))*self.CreationRate
			
			self.PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate/Norm
			self.PAnn = len(ExPairLocList)*self.AnnihilationRate/Norm
			self.PCre = (len(self.state) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*self.HoppingRate/Norm
			
			print RightHoppableLocList
			print LeftHoppableLocList
			print "**"
		
			r = self.NextActionProbability
			
			#print PHop, PAnn, PCre
			#print r
			#Hopping
			if r < self.PHop:
				HopSite = choice(RightHoppableLocList + LeftHoppableLocList)
				self.state[HopSite] = 0
				if HopSite in RightHoppableLocList:
					self.state[(HopSite+1)%len(self.state)] = 1
				else:
					self.state[(HopSite+1)%len(self.state)] = 0
					self.state[(HopSite)] = 1
					
				self.error_state[HopSite] = (self.error_state[HopSite] + 1) % 2
				#print "Hopping!"
				#print chain
		
			#Annihilating
			if (r >= self.PHop and r < self.PHop + self.PAnn):
				AnnihilateSite = choice(ExPairLocList)
				self.state[AnnihilateSite] = 0
				self.state[(AnnihilateSite+1)%len(self.state)] = 0
				
				self.error_state[AnnihilateSite] = (self.error_state[AnnihilateSite] + 1) % 2
				#print "Annihilating!"
				#print chain
		
			#Creating
			if (r >= self.PHop + self.PAnn):
				CreateSite = choice(EmptyPairLocList)
				self.state[CreateSite] = 1
				self.state[(CreateSite+1)%len(self.state)] = 1
				self.error_state[CreateSite] = (self.error_state[CreateSite] + 1) % 2
				#print "Creating!"
				#print chain
		
			ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(anyons_list)
			Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate + len(ExPairLocList)*self.AnnihilationRate + \
			(len(EmptyPairLocList))*self.CreationRate
			
			self.PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate/Norm
			self.PAnn = len(ExPairLocList)*self.AnnihilationRate/Norm
			self.PCre = (len(anyons_list) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*self.HoppingRate/Norm
			
			#Update the system time, next action time, and next action probability
			self.TotalTime = self.NextActionTime
			print "Action too late!" + str(self.TotalTime)
			self.NextActionProbability = random.random()
			self.NextActionTime = self.TotalTime + (-1./Norm)*np.log(self.NextActionProbability)
		else:
			#If we haven't exceeded the bath interaction timescale, we have to apply some swaps!
			anyons_list, CycleTime, NewRates, NoHops, Proceed, self.Sector = self.CorrectionProtocol(anyons_list, self.TotalTime, self.TotalTime+(1./self.CorrectionRate), self.CorrectionRate, \
																  self.PHop, self.PAnn, self.PCre, [action], self.Sector)
			#self.TotalTime+=CycleTime
			self.TotalTime+=1.
			print self.TotalTime
																  
			if NewRates == True:
				ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(anyons_list)
				Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate + len(ExPairLocList)*self.AnnihilationRate + \
				(len(EmptyPairLocList))*self.CreationRate
				
				self.PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*self.HoppingRate/Norm
				self.PAnn = len(ExPairLocList)*self.AnnihilationRate/Norm
				self.PCre = (len(anyons_list) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*self.HoppingRate/Norm
						
				self.NextActionProbability = random.random()
				self.NextActionTime = self.TotalTime + (-1./Norm)*np.log(self.NextActionProbability)
			
			self.state = anyons_list
		
		#Update the sector
		self.Sector = (self.Sector + self.CheckSector(self.state,p,pl,pr))%2
		
		
		
		done =  self.TotalTime > 1000. \
				or self.CheckState(self.state, self.Sector) == 1
		done = bool(done)

		if not done:
			reward = 1.0
		elif self.steps_beyond_done is None:
			# Pole just fell!
			self.steps_beyond_done = 0
			reward = 1.0
		else:
			if self.steps_beyond_done == 0:
				logger.warn("You are calling 'step()' even though this environment has already returned done = True. You should always call 'reset()' once you receive 'done = True' -- any further steps are undefined behavior.")
			self.steps_beyond_done += 1
			reward = 0.0

		return np.array(self.state), reward, done, {}

	def _reset(self):
		self.state = np.zeros(self.SystemLength)
		self.error_state = np.zeros(self.SystemLength)
		self.TotalTime = 0.
		self.state[0] = 1.
		self.state[22] = 1.
		self.error_state[22]=1
		self.error_state[23]=1
		self.steps_beyond_done = None
		return np.array(self.state)

	def _render(self, mode='human', close=False):
		if close:
			if self.viewer is not None:
				self.viewer.close()
				self.viewer = None
			return

		screen_width = 600
		screen_height = 400

		#world_width = self.x_threshold*2
		#scale = screen_width/world_width
		#carty = 100 # TOP OF CART
		polewidth = 10.0
		#polelen = scale * 1.0
		#cartwidth = 50.0
		#cartheight = 30.0
		

		if self.viewer is None:
			from gym.envs.classic_control import rendering
			self.viewer = rendering.Viewer(screen_width, screen_height)#, display=self.display)
			
			'''l,r,t,b = -cartwidth/2, cartwidth/2, cartheight/2, -cartheight/2
			axleoffset =cartheight/4.0
			cart = rendering.FilledPolygon([(l,b), (l,t), (r,t), (r,b)])
			self.carttrans = rendering.Transform()
			cart.add_attr(self.carttrans)
			self.viewer.add_geom(cart)
			l,r,t,b = -polewidth/2,polewidth/2,polelen-polewidth/2,-polewidth/2
			pole = rendering.FilledPolygon([(l,b), (l,t), (r,t), (r,b)])
			pole.set_color(.8,.6,.4)
			self.poletrans = rendering.Transform(translation=(0, axleoffset))
			pole.add_attr(self.poletrans)
			pole.add_attr(self.carttrans)
			self.viewer.add_geom(pole)'''
			
			for i in xrange(self.SystemLength):
			
				self.offsettrans = rendering.Transform()
				self.error_offsettrans = rendering.Transform()
				
				self.translist.append(self.offsettrans)
				self.error_translist.append(self.error_offsettrans)
				
				axle = rendering.make_circle(polewidth/2)
				error = rendering.make_circle(polewidth/4)
				
				axle.add_attr(self.offsettrans)
				error.add_attr(self.error_offsettrans)
				
				axle.set_color(.8,.6,.4)
				error.set_color(.1,.1,.1)
				
				self.viewer.add_geom(axle)
				self.viewer.add_geom(error)
				
				#print "Putting on the screen!"
			
			#self.track = rendering.Line((0,carty), (screen_width,carty))
			#self.track.set_color(0,0,0)
			#self.viewer.add_geom(self.track)
		for i,t in enumerate(self.translist):
			#print "something happening?"
			if self.state[i]!=0:
				#print "Moving to be visible!"
				t.set_translation(i*(400./self.SystemLength)+100., 200)
			else:
				t.set_translation(-10,-10)
		for i,t in enumerate(self.error_translist):
			if self.error_state[i]!=0:
				t.set_translation(i*(400./self.SystemLength)+100. + (400./self.SystemLength)/2., 200)
			else:
				t.set_translation(-10,-10)
		#print "This is being run, though!"
		#x = self.state
		#cartx = x[0]*scale+screen_width/2.0 # MIDDLE OF CART
		#self.carttrans.set_translation(cartx, carty)
		#self.poletrans.set_rotation(-x[2])

		return self.viewer.render()#return_rgb_array = mode=='rgb_array')
		
		
	#ISING CODE
	#*****************************************************
		
	def ReturnExcitationInformation(self, chain):
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



	def CalculateProbabilities(self, chain, CreationRate, AnnihilationRate, HoppingRate):
		ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(chain)
		Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate + len(ExPairLocList)*AnnihilationRate + \
		(len(EmptyPairLocList))*CreationRate
		
		
		PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate/Norm
		PAnn = len(ExPairLocList)*AnnihilationRate/Norm
		PCre = (len(chain) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*HoppingRate/Norm
		
		return PHop, PAnn, PCre



	def AdvanceTime(self, chain, StartTime, CreationRate, AnnihilationRate, HoppingRate, sector):
		ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(chain)
		Norm = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate + len(ExPairLocList)*AnnihilationRate + \
		(len(EmptyPairLocList))*CreationRate
		
		PHop = (len(RightHoppableLocList)+len(LeftHoppableLocList))*HoppingRate/Norm
		PAnn = len(ExPairLocList)*AnnihilationRate/Norm
		PCre = (len(chain) - len(ExPairLocList) - (len(RightHoppableLocList)+len(LeftHoppableLocList)))*HoppingRate/Norm
				
		r = random.random()
		DeltaTau = (-1./Norm)*np.log(r)
		
		chain, CycleTime, NewRates, NoHops, Proceed, sector = CorrectionProtocol(chain, StartTime, StartTime+DeltaTau, CorrectionRate, \
															  PHop, PAnn, PCre, CorrectionSwaps, sector)
		#NewRates = False
		#CycleTime = 0
		#NoHops = True
		
		p = int(floor(len(chain)/2.))
		pl,pr = chain[p],chain[p+1] #previous values of chain
		
		if NewRates == False:
			ExLocList, ExPairLocList, EmptyLocList, EmptyPairLocList, RightHoppableLocList, LeftHoppableLocList = self.ReturnExcitationInformation(chain)
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
				chain[HopSite] = 0
				if HopSite in RightHoppableLocList:
					chain[(HopSite+1)%len(chain)] = 1
				else:
					chain[(HopSite-1)%len(chain)] = 1
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

	def CheckSector(self, chain,p,pl,pr):
		increment = 0
		
		if chain[p]!=pl and chain[p+1] != pr:
			increment = 1
		
		
		#print p,pl,pr,"\t",chain[p],chain[p+1],"\t",increment
		#print chain
		
		return increment


	#Constructs a list with the indices for conditional swaps in the correction protocol
	#Convention is that the value at protocol[i] is CSWAPPED with protocol[(i+1)%length]
	def SwapProtocol(self, length):
		sublength = length/2 - 1
		protocol = []
		for i in xrange(length):
			for j in xrange(sublength):
				for k in xrange(sublength - j):
					protocol.append((i+(j+k))%length)
				for k in xrange(sublength - j):
					protocol.append((i+(sublength-k-1))%length)
					
		return protocol



	def SwapProtocol2(self, length):
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

	def SwapProtocol3(self, length):
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

	def SwapProtocol4(self, length):
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

	def SwapProtocol5(self, length):
		sublength = int(math.ceil(length/2))
		protocol = []

		for j in xrange(sublength):
			if j%2==0:
				protocol.append(2*j)
				protocol.append(2*j+1)
				protocol.append(2*j)
				protocol.append(2*j+1)


		return protocol



	def CSwap(self, chain, i):
		#print i
		if chain[i]!=chain[(i+1)%len(chain)]:
			inter = chain[i]
			chain[i] = chain[(i+1)%len(chain)]
			chain[(i+1)%len(chain)] = inter
			self.error_state[i] = (self.error_state[i] + 1) % 2
		####print "Swapping at " + str(i) + "!: ",chain
		return chain
		
	def CorrectionProtocol(self, chain, oldtime, newtime, CorrectionRate, PHop, PAnn, PCre, CorrectionSwaps, sector):
		
		print "Attempting correction protocol"
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
		
		'''if (oldtime + CorrectionPeriod/NumberOfSwaps) > newtime:
			Proceed = False
			#psuccess = (newtime - oldtime) / (CorrectionPeriod/NumberOfSwaps)
			#if random.random() < psuccess:
			#	Proceed = True
		'''		
		#This loop should only execute once. lol why is it here then.  Because I might generalize this later
		print oldtime+CycleTime < newtime
		print not(ProbabilityHasChanged)
		print not(PHop == 0)
		while(oldtime+CycleTime < newtime and not(ProbabilityHasChanged) and not(PHop == 0)):# and Proceed == True):# and not(PAnn > 0)):

			####print "Timing information: ", CycleTime,"\t",oldtime,"\t", newtime,"\t",(newtime-oldtime)-CycleTime,"\t",CorrectionPeriod/NumberOfSwaps
			#chain = self.CSwap(chain, CorrectionSwaps[IndexInCycle])
			chain = self.CSwap(chain, CorrectionSwaps[0])
			print "Swapping" + str(CorrectionSwaps[0])
			#parallel?
			#chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + NumberOfSwaps/4)%NumberOfSwaps])
			#chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 2*NumberOfSwaps/4)%NumberOfSwaps])
			#chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 3*NumberOfSwaps/4)%NumberOfSwaps])
			
			#chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 2*NumberOfSwaps/4)%NumberOfSwaps])
			#chain = CSwap(chain, CorrectionSwaps[(IndexInCycle + 3*NumberOfSwaps/4)%NumberOfSwaps])
			
			PHopInter, PAnnInter, PCreInter = self.CalculateProbabilities(chain, self.CreationRate, self.AnnihilationRate, self.HoppingRate)
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

	def CheckState(self, chain, sector):
		if sum(chain)==0:
			return 2*sector-1
		else:
			return 0


	def ProcessTraj(self, traj,avgtraj,maxtime):
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
