"""Constants to be used as parameters for simulation"""
import numpy as np 

# Molecular Size
Size = {"N2": 2,
	"H2O": 3}


# Atomic Masses of Elements [amu]
AMU = {"H": 1.00784,
	"He": 4.002602,
	"Li": 6.938,
	"C": 12.0096,
	"N": 14.0067,
	"O": 15.999,
	"F": 18.998403,
	"S": 32.066,
	"Ne": 20.1797 }

Mass = {"N2": 2*AMU["N"],
	"H2O": AMU["O"]+2*AMU["H"]}


Initial_Position = {"H2O":np.array([ [0, -0.757, 0.757],     # current in the form [X, Y, Z]^T
				   [0.06556811, -0.52043189, -0.52043189],
				   [0, 0, 0]]),
		"N2": np.array([[0,0],[0.54875,-0.54875],[0,0]])}

Collision_Sphere = {"H2O": 1,
		"N2": 0.548}
																				
mol_type = {"H2O":["O", "H", "H"],
	    "N2": ["N", "N"]}

# TIP3P Data
TIP3P = {"OH_distance": 0.9572, # From TIP3P model
		 "HH_distance": 1.514, #Calculated using trig
		 "Dihedral": 104.52, # Dihedral angle
		 "A": 582.0,		 # Parameter A for Lennard Jones
		 "B": 595.0,		 # Parameter A for Lennard Jones
		 "q_O": -0.834,		 # Charge of oxygen
		 "q_H": 0.417,		 # Charge of hydrogen
		 "k_c": 332.1        #[A kcal/ mol e^2] Electrostatic constant
}

# System Temperature
T = 1.5

# Maximum step size for normal steps (Angstroms)
L = 1.0

# Magnification factor for L
beta = 3.0

# Probability of a translational magstep
P_MT = 0

# Probability of a rotational magstep
P_RT = 0

# Number of Monte Carlo moves
M = 100

# Molecule Size: 3 for triatomics, 2 for diatomics
mol_size = 2

# Number of Molecules
N = 50

k = 1.38064853e-23 * 6.022140858e23 * 1e-3  # [kJ/mol-K]  Boltzmann constant (J/K converted using Avogadro constant, and 1/1000 for kJ)
sigma = 3.418
SIGSQ = sigma**2
epsilon = 124*k # [kJ/mol] From McQuarrie pg. 434

L_MAG = beta * L

# Import data from a text file. Cancelling this to place a widget instead.
file = open("./user_inputs.txt")
f = file.readlines()

# N = int(f[1])
mol_size = int(f[2])

file.close()

