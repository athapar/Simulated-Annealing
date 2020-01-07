""" Functions to transform molecules while Mag-Walking.
	
	Includes linear translation of molecules and rotations 

	Array Shape:

	    1	2	3  ...  n
	x:[x1  x2   x3 ...  xn]
	y:[y1  y2   y3 ...  yn]
	z:[z1  z2   z3 ...  zn]


	Particle positions in array:

	x_i == R[i][0]
	y_i == R[i][1] 
	z_i == R[i][2]"""

import numpy as np
import math
import parameters as pm


# T = 1.5

# # Maximum step size for normal steps (Angstroms)
# L = 2.0

# # Magnification factor for L
# beta = 3.0

# # Probability of a translational magstep
# P_MT = 0.2

# # Probability of a rotational magstep
# P_RT = 0.3

# # Number of Monte Carlo moves
# M = 100

# # Molecule Size: 3 for triatomics, 2 for diatomics
# mol_size = 2

# # Number of Molecules
# N = 2

# k = 1.38064853e-23 * 6.022140858e23 * 1e-3  # [kJ/mol-K]  Boltzmann constant (J/K converted using Avogadro constant, and 1/1000 for kJ)
# sigma = 3.418
# SIGSQ = sigma**2
# epsilon = 124*k # [kJ/mol] From McQuarrie pg. 434

# L_MAG = beta * L




# Potential Calculation
def LJ_Potential(size, num_particles, X,Y,Z):
	V = 0.0
	print(X)
	for a in range(size):
		for b in range(size):
			for i in range(num_particles-1):
				for j in range(i+1,num_particles):
					RXAB = X[i,a] - X[j,b]
					RYAB = Y[i,a] - Y[j,b]
					RZAB = Z[i,a] - Z[j,b]

					RIJSQ = RXAB**2 + RYAB**2 + RZAB**2
					SR2 = pm.SIGSQ/RIJSQ
					SR6 = SR2**3
					SR12 = SR6**2
					V += SR12 - SR6



	V *= 4*pm.epsilon
	return V

# Check the rotation functions again since the molecules aren't changing their X-axis positions
# So probably not the xrotation function as it does not change x positions, and since y and z positions do change
# Z positions are also looking the same


# Object oriented version
# Takes inputs in degrees (will convert to radians inside the function)
def xrotation(theta, pos, com_frame):
	x = pos
	theta = theta*math.pi/180 # converts theta to radians
	xvec = np.array([[1,0,0],
			[0,math.cos(theta),math.sin(theta)],
			[0,-math.sin(theta), math.cos(theta)]])

	R = np.dot(xvec,com_frame) # Find new matrix in COM Frame
	Delta = R - com_frame 	   # Find changes in position from rotation

	x += Delta      			# Add difference to position matrix
	return(pos)


def yrotation(theta, pos, com_frame):
	theta = theta*math.pi/180  # Converts angle to radians
	yvec = np.array([[0,math.cos(theta),math.sin(theta)],
				 	 [0,1,0],
			     	 [-math.sin(theta),0, math.cos(theta)]])

	R = np.dot(yvec,com_frame) # Find new matrix in COM Frame
	Delta = R - com_frame 	   # Find changes in position from rotation

	pos += Delta      			# Add difference to position matrix
	return(pos)


def zrotation(theta, pos, com_frame):
	theta = theta*math.pi/180  # Converts theta to radians
	zvec = np.vstack([[math.cos(theta), -math.sin(theta), 0],
					 [math.sin(theta),math.cos(theta),0],
					 [0,0,1]])

	R = np.dot(yvec,com_frame) # Find new matrix in COM Frame
	Delta = R - com_frame 	   # Find changes in position from rotation

	pos += Delta      			# Add difference to position matrix
	return(pos)

def translation(pos, delta):
	pos += np.array([delta[0]*np.ones(pos.shape[0]),
					delta[1]*np.ones(pos.shape[0]),
					delta[2]*np.ones(pos.shape[0])
])
	return pos 


# Translations
def xtranslation(pos, delx):
	pos[0] += delx*np.ones(pos.shape[0])
	return pos
def ytranslation(pos, dely):
	pos[1] += dely*np.ones(pos.shape[0])
	return pos
def ztranslation(pos, delz):
	pos[2] += delz*np.ones(pos.shape[0])
	return pos
