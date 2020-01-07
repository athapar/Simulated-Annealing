""" Input Array Shape:
    1	2	3  ...  n
x:[x1  x2   x3 ...  xn]
y:[y1  y2   y3 ...  yn]
z:[z1  z2   z3 ...  zn]


Particle positions in array:
x_i == R[i][0]
y_i == R[i][1] 
z_i == R[i][2]

One issue to look out is that altering copy may affect original variable. Using deepcopy works but don't know about speed. X.copy() and deepcopy(X) both work
https://stackoverflow.com/questions/8771808/copy-a-dictionary-into-a-new-variable-without-maintaining-the-link-with-previous

There is an issue with y-rotations and multiple rotation of a molecule.
Perform hand calculation to check
Main issue is that rigid molecule is not maintained. OH distance is changed
"""




## %% Import modules and define parameters
import numpy as np
import parameters as pm
import math
import Calculations as calc
import data_handler as data
import re

# Basic idea is to make each particle have a matrix of its own positions and create a simple loop to add them to
# the x,y,z matrices for the system. Will need to check molecules aren't overlapping at the beginning of the loop
# First construct setup with some basic properties: Size, Position Matrix, Can add tags for printing (more complicated so later) 

# Want a widget system to select molecules

# First create an initial position for the molecule with COM at 0, make it the default position
# Then perform random operations for rotation and translation to move it in space
# Want to move COM and maintain intramolecular distances when this is done
# Lastly, add a checking system to make sure molecules are not overlapping
# This will be difficult because how do you know molecules overlap in 3d?


MAXT = 5 # Maximum translation
MAXR = 360.0 # Maximum rotation
COL_DIST = 2 # [A] Collision distance
MAX_SEARCH = 100 # Max number of attempts to find a position that isn't colliding with other particles
INCREASE_MOVE = 1.3 # Factor to increase move size if algorithm fails to find position without colliding
NAME = "H2SO4"

# Creates a molecule with initial position given by self.position under the __init__() function
# Each molecule has functions to find COM and positions with COM as origin(COM_FRAME()) for rotations


## %% Define individual molecules
class molecule():
	# Initial Parameters
	def __init__(self, name, custom_path="Sim/Custom.xyz"):
		self.name = name # Specify type of molecule by name
		if self.name.upper() == "H2O" or self.name.upper() == "WATER":
			# Basic parameters
			self.name = "H2O"
			self.size = 3 						# Number of atoms
			self.mol_type = ["O", "H", "H"] 	    # Keys for atoms
			self.col_sphere = 2 				# Collision-sphere [A]
			self.mass = 2*pm.AMU["H"] + pm.AMU["O"] # Molecular mass
			self.center = [0,0,0] # Initial center is origin

			self.position = np.array([ [0, -0.757, 0.757],     # current in the form [A1 A2 A3...] where Ai = [xi yi zi]^T
					   [0.06556811, -0.52043189, -0.52043189],
					   [0, 0, 0] ])
			
			# Make one copy the COM Frame of the molecule
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()


		elif self.name.upper() == "N2" or self.name.upper() == "NITROGEN":
			self.name = "N2"
			# Basic Parameters
			self.size = 2  			# Number of atoms
			self.mol_type = ["N", "N"]	# Keys for atoms
			self.col_sphere = 0.54875        # Collision-sphere [A]
			self.mass = 2*pm.AMU["N"]
			self.position = np.array([[0,0],[0.54875,-0.54875],[0,0]]) # Initialize position
			self.center = [0,0,0]
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()
			
		elif self.name.upper() == "NH3" or self.name == "ammonia" or self.name == "Ammonia":
			self.name = "NH3"
			
			# Basic Parameters
			self.size = 4
			self. mol_type = ["N", "H", "H", "H"]
			self.col_sphere = 3 
			self.mass = pm.AMU["N"] + 3*pm.AMU["H"]
			
			self.position = np.array([ [0, 0, 0.8068, -0.8068],
							 [0, 0.9316, -0.4658, -0.4658],
							 [0.1111, -0.2592, -0.2592, -0.2592] ])
	
			self.center = [0,0,0.04535749]
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()
			
		elif self.name.upper() == "SF4":
			self.size = 5
			# Basic Parameters
			self.mol_type = ["S", "F", "F", "F", "F"]
			self.col_sphere = 4
			self.mass = pm.AMU["S"] + 4*pm.AMU["F"]
			self.position = np.array([ [0, 0, 0, 1.2055, -1.2055],
								       [0, 1.6255, -1.6255, 0, 0],
								       [0.3835, 0.2401, 0.2401, -0.5810, -0.5810] ])
			self.center = [0,0,-0.00606887]
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()
		
		elif self.name.upper() == "SF6":
			self.size = 7
			# Basic Parameters
			self.mol_type = ["S", "F", "F", "F", "F", "F", "F"]
			self.col_sphere = 6
			self.mass = pm.AMU["S"] + 6*pm.AMU["F"]
			self.position = np.array([ [0, 0, 0, 1.554, 0, -1.554, 0],
								       [0, 0, 1.554, 0, -1.554, 0, 0],
									    [0, 1.554, 0, 0, 0, 0, -1.554],])
			self.center = [0,0,0]
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()
			
		elif self.name.upper() == "H2SO4":
			self.size = 7
			# Basic Parameters
			self.mol_type = ["S", "O", "O", "O", "O", "H", "H"]
			self.col_sphere = 6
			self.mass = pm.AMU["S"] + 4*pm.AMU["O"] + 2*pm.AMU["H"]
			self.position = np.array([[0, 0, 0, 1.2193	, -1.2193, 1.4709, -1.4709],
									 [0, 1.2438, -1.2438, 0.0242, -0.0242, -0.8641, 0.8641],
									 [0.1534, 0.8192	, 0.8192	, -0.8376, -0.8376, -1.0800, -1.0800	]])	
			self.center = [0,0,0]
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()
		
		# 
		elif self.name.upper() == "CUSTOM":
			with open(custom_path, "r") as f:
				f1 = f.readlines()
				self.size = int(f1[0].split()[0])
				self.name = f1[1].split()[0]
				
				self.mol_type = []
				self.mass = 0
	
				for i in range(2, self.size+2):
					Split = f1[i].split()
					atom_name = Split[0]
					self.mol_type.append(atom_name)
					self.mass += pm.AMU[atom_name]
					if i == 2:
						self.position = np.array([float(i) for i in Split[1:]])
					else:
						self.position = np.vstack([self.position,np.array([float(i) for i in Split[1:]])])

#			self.position = self.position.transpose()
			self.center = [0,0,0]
			self.com_frame = self.position.copy()
			self.previous_position = self.position.copy()
			self.col_sphere = 6
			

	def COM(self):
		R = self.position.copy()
		self.center = [0,0,0]
		A = 0
		for i in range(self.size):
			A += pm.AMU[self.mol_type[i]] * np.array([j for j in R[:,i]]) # m_i * r_i (Use AMU to get m_i and self.postion list comprehension to find r_i)
		self.center = A/self.mass
		return self.center

	def COM_FRAME(self):
		self.com_frame = 0
		R = self.position.copy()
		for i in range(3):
#			print(x[i])
			R[i] -= self.center[i]*np.ones(self.size)
#			print(x[i])
		self.com_frame = R
		return self.com_frame


	def translate(self,delta):
		
		self.position += np.array([delta[0]*np.ones(self.position.shape[1]),
					delta[1]*np.ones(self.position.shape[1]),
					delta[2]*np.ones(self.position.shape[1])
					])

		# Update COM and COM
		return self.position

	def xrotation(self, theta):
		theta = theta*math.pi/180 # converts theta to radians
		xvec = np.array([[1,0,0],
				[0,np.cos(theta),np.sin(theta)],
				[0,-np.sin(theta), np.cos(theta)]])

		R = np.dot(xvec,self.com_frame) # Find new matrix in COM Frame
		Delta = R - self.com_frame	   # Find changes in position from rotation

		self.position += Delta      			# Add difference to position matrix
		return(self.position)

	def yrotation(self, theta):
		theta = theta*math.pi/180 # converts theta to radians
		yvec = np.array([[np.cos(theta), 0 , np.sin(theta)],
					 	 [0,1,0],
				     	 [-np.sin(theta),0, np.cos(theta)]])

		R = np.dot(yvec,self.com_frame) # Find new matrix in COM Frame
		Delta = R - self.com_frame	   # Find changes in position from rotation

		self.position += Delta      			# Add difference to position matrix
		return(self.position)

	def zrotation(self, theta):
		theta = theta*math.pi/180 # converts theta to radians
		zvec = np.vstack([[np.cos(theta), -np.sin(theta), 0],
							 [np.sin(theta),np.cos(theta),0],
							 [0,0,1]])

		R = np.dot(zvec,self.com_frame) # Find new matrix in COM Frame
		Delta = R - self.com_frame	   # Find changes in position from rotation

		self.position += Delta      			# Add difference to position matrix

		return(self.position)
		
	def update(self):
		self.center = self.COM()
		self.com_frame = self.COM_FRAME()		
		return self.position

	# Copies previous position and then randomly moves the particle
	def random_move(self, dist=1, theta=90):
		self.xrotation(theta*np.random.random())
		self.update()
		
		self.yrotation(theta*np.random.random())
		self.update()
		
		self.zrotation(theta*np.random.random())
		self.update()
		
		self.translate([dist*(2*np.random.random()-1),  # Translate each randomly
				dist*(2*np.random.random()-1),
				dist*(2*np.random.random()-1)]) 
		self.update()
		
		return self.position



## %% Use molecule class to create initial configuration of the system.
		
def create_system(Npart, Mol_Class=molecule, name=NAME, create=True, filename="testprint"):
	M=[]
	for i in range(Npart):
		M.append(Mol_Class(name))  # Initialize the molecule	
		M[i].random_move(MAXT, MAXR) # Randomly move particle

		# This portion ensures particle 
		if i > 0: # Start checking after 1st particle
			j = 0 # Initialize j for while loop
			max_search = 100
			move_size = MAXT # Initial move_size is set by MAXT, can later be increased
			m = 0
			while j < i:
				if calc.molecular_dist(M[i],M[j]) >= 2*COL_DIST: # Check if within collision distance
					j += 1
				else:
					M[i].position = M[i].previous_position # Set equal to previous position
					M[i].update() # Update COM and coordinates with COM at origin
					M[i].random_move(move_size, MAXR) # Randomly move again
					j = 0
					m += 1
				if m > max_search:
					move_size = INCREASE_MOVE*move_size
					m = 0
	if create:
		data.to_xyz(M, filename)
	return M


#
#def custom(file="Sim/Custom.xyz"):
#	f = open(file, "r")
#	f1 = f.readlines()
#	Size = int(f1[0].strip())
#	Name = f1[1].strip()
#	
#	Split = re.split(r"\t",Strip)
#	name = Split[0].strip()
#	mol_type = [name]
#	mass = pm.AMU[name]
#	position = np.array([float(i) for i in Split[1:]])
#	
#	for i in range(2, len(f1)):
#		Strip = f1[i].rstrip()
#		Split = re.split(r"\t",Strip)
#		name = Split[0].strip()
#		mol_type.append(name)
#		mass += pm.AMU[name]
#		position = np.vstack([position,np.array([float(i) for i in Split[1:]])])
#	
#	return Size, Name, position, mass, mol_type
#		
	
	
#Size, Name, position, mass, mol_type = custom()


#distance = molecular_dist(sys[0],sys[1])
#df = create_df(sys)

