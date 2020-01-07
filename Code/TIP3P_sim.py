"""
TIP3P Simulation
Generates and perturbs system of water molecules according to TIP3P paramters. 
"""
import numpy as np
import pandas as pd
import Initial_Configuration as IC
import parameters as pm
from tqdm import tqdm
import data_handler as data
import Calculations as calc

# Meeting 10-27
# Don't record every move in the algorithm
# Check the acceptance ratio and make sure it's not too high
# Change the structure of the data so that it is [X Y Z] where each is a column vector


A = pm.TIP3P["A"]
B = pm.TIP3P["B"]
k_c = pm.TIP3P["k_c"]

CUT_OFF = 50 # [A] cut-off distance for interactions
file = "TIP3P"    # File name to be saved in Sim directory
N_part = 10        # Number of particles
Col_Dist = 1 #[A]
k_B = 0.0019872041 # [kcal/mol-K]

T = 10            #[K]

NAME = "H2O"
THETA = 20
SAVE = 90

class Water_Mol(IC.molecule):
	# Change initialization to include give charges and specify oxygen and hydrogen positions separately
	def __init__(self, name = "H2O"): 
		super().__init__("H2O") # Gives default values for water
		
		self.name = name
		# Add TIP3P parameters
		# Add charges for coulombic force for hydrogen and oxygen 
		self.charges = [pm.TIP3P["q_O"], pm.TIP3P["q_H"], pm.TIP3P["q_H"]] # Arrangement is [O H H]

		# Keep track of the atom positions 
		self.Ox = self.position[:,0]
		self.H1 = self.position[:,1]
		self.H2 = self.position[:,2]
		# Test position for particle
		self.previous_position = self.position.copy()
		
	def update(self):
		super().update()
		
		# Change the update function to include updating Oxygen and Hydrogen positions
		self.Ox = self.position[:,0]
		self.H1 = self.position[:,1]
		self.H2 = self.position[:,2]	
		
	def random_rotation(self, angle):
		self.xrotation(angle)
		self.update()
		
		self.yrotation(angle)
		self.update()
		
		self.zrotation(angle)
		self.update()
		
		return self.position
		
		
		
	
# Use potential for 3-site model
def TIP3P_Potential(M):
	E = 0
	for i in range(len(M)-1):
		for j in range(i+1,len(M)):
			r_oo = calc.distance(M[i].Ox, M[j].Ox)
			EIJ = 0
			r6 = r_oo ** 6
			r12 = r6 ** 2
			for a in range(3):
				for b in range(3):
					EIJ += M[i].charges[a] * M[j].charges[b] / calc.distance(M[i].position[:,a],M[j].position[:,b])
			EIJ = k_c * EIJ
			EIJ = EIJ + A/r12 - B/r6

		for k in range(3):
			E += (M[i].center[k])**2
	return E



def Simulation(N_part,trials, T=10, n_eq = 1000, save=10):
	# Initialize arrays
	M_array = [] # Array stores all snapshots
	E_array = []
	E2_array = []
	
	# Create random initial position
	M = IC.create_system(N_part, name="Water", Mol_Class=Water_Mol, filename= file)
	X, Y, Z = data.Store_Data(M)
	E = TIP3P_Potential(M)
	
	# Create XYZ File
	data.to_xyz(M, filename = file)
	
	M_array.append(M)
	E_array.append(E)
	E2_array.append(E**2)
	
	accepted = 0
	rejected = 0
	total = 0
	

	for k in tqdm(range(trials)):
		for i in range(N_part):
			
			M[i].previous_position = M[i].position
			
			# Translation condition for MagWalk
			if np.random.random() < pm.P_MT:
				vector = pm.L_MAG * np.array([ (2*np.random.random()-1), (2*np.random.random()-1), (2*np.random.random()-1)  ])
				M[i].translate(vector)
			else:
				vector = pm.L * np.array([ (2*np.random.random()-1), (2*np.random.random()-1), (2*np.random.random()-1)  ])
				M[i].translate(vector)
			
			# Rotation condition for MagWalk
			if np.random.random() < pm.P_RT:
				M[i].random_rotation(180)
			else:
				M[i].random_rotation(THETA)

			
			Collide = False
			# Loop to check if molecule collides with other molecules 
			for j in range(N_part):
				if j == i:
					pass
				elif calc.distance(M[i].Ox, M[j].Ox) < 2*Col_Dist:
					Collide = True
					break
				else:
					pass
			
			# REJECT MOVE BASED ON COLLISION
			if Collide:

				if  k> n_eq:
					rejected += 1
					total += 1
					
					testE = TIP3P_Potential(M)	
					
					E_array.append(testE)
					E2_array.append(testE**2)
					M_array.append(M)
				
				X, Y, Z = data.Append_Data(M, X, Y, Z)
				M[i].position = M[i].previous_position
				print("REJECT: Collide ") 
				
				if k%save == 0:
					data.append_xyz(M, filename = file)
					
				
				
				
			# IF NOT COLLIDING, BOLTZMANN TEST
			else:
				testE = TIP3P_Potential(M)
				deltaE = testE - E
				
				# ENERGY DECREASE = ACCEPT
				if np.sign(deltaE) == -1 or np.sign(deltaE) == 0:
					
					X, Y, Z = data.Append_Data(M, X, Y, Z)
					print("ACCEPT: ENE-D \t E=",E,"test=", testE)
					
					
					if k%save == 0:
						data.append_xyz(M, filename = file)
						
					if k > n_eq:
						accepted += 1
						total += 1
						
						E_array.append(testE)
						E2_array.append(testE**2)
						M_array.append(M)
	#					
					E = testE
					
				# PASSES BOLTZMANN TEST = ACCEPT
				elif np.random.random() <= np.exp(-deltaE/k_B/T):
						
					if k > n_eq:
						accepted += 1
						total += 1
						
						E_array.append(testE)
						E2_array.append(testE**2)
						M_array.append(M)
						
					E = testE
					X, Y, Z = data.Append_Data(M, X, Y, Z)
					print("ACCEPT: PASSB \t E=",E,"test=", testE, "Bolz", np.exp(-deltaE/k_B/T))
					
					if k%save == 0:
						data.append_xyz(M, filename = file)
						
						
				# FAILS BOLTMANN TEST = REJECT
				else:
						
					if k > n_eq:
						rejected += 1
						total += 1
	
						E_array.append(testE)
						E2_array.append(testE**2)
						M_array.append(M)

					X, Y, Z = data.Append_Data(M, X, Y, Z)
					data.append_xyz(M, filename = file)
					M[i].position = M[i].previous_position

					if k%save == 0:
						print("REJECT: FAILB \t E=",E,"test=", testE, "Bolz", np.exp(-deltaE/k_B/T))					
				
				
				
			# WANT TO NOT CONTINUE IF PARTICLE COLLIDES - LOOK AT OLD CODE
			# Will try to use Boolean logic and see how it works
		
				
			# CHECK IF PARTICLES COLLIDE - Can immediately reject
			# CHECK POTENTIAL CHANGE
			
			# IF WORKS store data
			# IF NOT M[i].position = M[i].previous_position
			
	# Change indices to start from 1 instead of 0 (python default)
	X.index = range(1, len(X.index)+1)
	X.columns = range(1, len(X.columns)+1)
	
	Y.index = range(1, len(Y.index)+1)
	Y.columns = range(1, len(Y.columns)+1)
	
	Z.index = range(1, len(Z.index)+1)
	Z.columns = range(1, len(Z.columns)+1)
	
	
	return X, Y, Z, np.array(E_array), np.array(E2_array), np.array(M_array), accepted, rejected, total
	

X, Y, Z, E_array, E2_array, M_array, accepted, rejected, total = Simulation(2, 100, n_eq=0, save=1)
print("Acceptance:", accepted/total*100,"%")

#X, Y, Z, E_array, M_array = np.zeros(10), np.zeros(10), np.zeros(10), np.zeros(10), np.zeros(10)
#
#for i in range(10):
#	X[i], Y[i], Z[i], E_array[i], E2_array[i], M_array[i] = Simulation(10,20, T=i)
#T1=1
#T2=2
#T3=3
#
#X1, Y1, Z1, E_array1, E2_array1, M_array1, accepted1, rejected1, total1 = Simulation(5, 10000, n_eq=3000, T=T1)
#X2, Y2, Z2, E_array2, E2_array2, M_array2, accepted2, rejected2, total2 = Simulation(5, 10000, n_eq=3000, T=T2)
#X3, Y3, Z3, E_array3, E2_array3, M_array3, accepted3, rejected3, total3 = Simulation(5, 10000, n_eq=3000, T=T3)



#import matplotlib.pyplot as plt
#
#plt.plot(E_array3[1:])
#plt.ylabel("Energy")
#plt.title("Energy vs Trials")
#plt.show()

# Initialize System
def CREATE_XYZ(MM,file):
	for i in range(len(MM)):
		if i == 0:
#			to_xyz()
			pass
		else:
#			append_xyz()
			pass

			
		

