"""Constants to be used as parameters for simulation"""
import numpy as np 
import csv

# Physical Constants
K_C = 332.0626
K_B = 0.0019872041 # Boltzmann constant [kcal/mol-K]
R = 1.987e-3 # Gas constant [kcal/mol/K]
HARTREE_TO_KCAL = 627.50908 # Kcals per Hartree
BOHR_TO_ANGSTROM = 0.529177 # Angstroms per Bohr radius


# Default Inputs for functions
SIGMA = 3.0 
EPSILON = 0.2 


# Ammonium Halide Sim Data
# Taken from  Topper 2019 OPLS parameters file in OneDrive
# Data for JB/RQT Potential and 1996 j-walking paper
AmmCl = {"A"  :{"N":{"N":65725.30104, "H":6474.638687, "Cl":71674.08712},
                "H": {"N":6474.638687, "H":637.6747271, "Cl":0},
                "Cl": {"N":71674.08712, "H":0, "Cl":78156.25591}
                },
         
         "B": {"N":{"N":2.950052629, "H":3.359934389, "Cl":3.127497983},
               "H": {"N":3.359934389, "H":3.770005121, "Cl":0},
               "Cl": {"N":3.127497983, "H":0, "Cl":3.304943337}
               },
         
         "C": {"N":{"N":349.8966515, "H":120.1950735, "Cl":740.4421086},
               "H": {"N":120.1950735, "H":41.30056447, "Cl":138.247277},
               "Cl": {"N":740.4421086, "H":138.247277,"Cl":1566.425839}
               },
         "D": {"N":{"N":15.12869743, "H":19.66730666, "Cl":60.51478972},
               "H": {"N":19.66730666, "H":611.5049759, "Cl":13278.15516},
               "Cl": {"N":60.51478972, "H":13278.15516, "Cl":242.0591589}
               },
         "Q": {"N":{"N":0.16, "H":-0.14, "Cl":0.4},
               "H": {"N":-0.14, "H":0.1225, "Cl":-0.35},
               "Cl": {"N":0.4, "H":-0.35, "Cl":1}
               } # Do not need comma for last element in dictionary, but do to separate other elements
            }


# Collision spheres currently used for molecules
Collision_Sphere = {"H2O": 2,
		"N2": 0.548,
		"SO2": 2.47,
		"NH3": 3.0,
		"SF4": 4.0,
		"SF6": 3.1214+1,
		"H2SO4": 4.0492+1,
		"NH4": 1.68,
        "Cl": 1.75,
        "HCl": 1.2746+1.02,
        "F": 1.47,
        "Br": 1.85,
        "I": 1.98
		}
																				
# TIP3P Simulation Data 
TIP3P = {"OH_DISTANCE": 0.9572, # From TIP3P model
		 "HH_DISTANCE": 1.514, # Calculated manually OH distance and dihedral
		 "DIHEDRAL": 104.52, # Dihedral angle
		 "A": 582.0*1000,		 # [kcal A^12/mol]Parameter A for Lennard Jones
		 "B": 595.0,		 #[kcal A^6/mol] Parameter B for Lennard Jones
		 "SIGMA":3.1506,
		 "EPSILON":76.54*1.987e-3,
		 "Q_O": -0.834,		 # Charge of oxygen
		 "Q_H": 0.417,		 # Charge of hydrogen
		 "K_C": 332.1        #[A kcal/ mol e^2] Electrostatic constant
}

# TIP4P Simulation Data
TIP4P = {"r_OH": 0.9572,
		 "r_OM": 0.15,
		 "HOH_angle": 104.52,
		 "sigma": 3.154,
		 "A": 600.0e-3,
		 "B": 610.0,
		 "q_M": -1.04,
		 "q_H": 0.52,
		 "k_c": 332.1        #[A kcal/ mol e^2] Electrostatic constant
		 }


# sigma = {"S": {"S":3.39, "O":41.195},
#          "O": {"S":41.195, "O":79}}

# epsilon = {"S": {"S":73.8, "O":15.0029997},
#          "O": {"S":15.0029997, "O":3.05}}

# Lennard Jones parameters for other molecules
SO2 = {"S": {"epsilon": 73.8, "sigma": 3.39, "q": 0.59},
       "O": {"epsilon": 79.0, "sigma": 3.05, "q": -0.295},
       
       "sigma": {"S": {"S":3.39, "O":41.195},
                 "O": {"S":41.195, "O":79}},
       
       "epsilon": {"S": {"S":73.8, "O":15.0029997},
                   "O": {"S":15.0029997, "O":3.05}}
       }



NH3 = {"N": {"epsilon": 85.458, "sigma": 3.42, "q": -1.02},
	     "H": { "q": 0.34}}

CO2 = { "C": {"epsilon": 52.84, "sigma": 3.75, "q": 0.7},
	     "O": {"epsilon": 63.41, "sigma": 2.96, "q": -0.35}}



# Atomic Lennard Jones Parameters
# Table E1. https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527676750.app5


# Atomic Masses of Elements [amu]
AMU = {'H': 1.00794,
		 'He': 4.002602,
		 'Li': 6.941,
		 'Be': 9.01218,
		 'B': 10.811,
		 'C': 12.011,
		 'N': 14.00674,
		 'O': 15.9994,
		 'F': 18.998403,
		 'Ne': 20.1797,
		 'Na': 22.989768,
		 'Mg': 24.305,
		 'Al': 26.981539,
		 'Si': 28.0855,
		 'P': 30.973762,
		 'S': 32.066,
		 'Cl': 35.4527,
		 'Ar': 39.948,
		 'K': 39.0983,
		 'Ca': 40.078,
		 'Sc': 44.95591,
		 'Ti': 47.88,
		 'V': 50.9415,
		 'Cr': 51.9961,
		 'Mn': 54.93805,
		 'Fe': 55.847,
		 'Co': 58.9332,
		 'Ni': 58.6934,
		 'Cu': 63.546,
		 'Zn': 65.39,
		 'Ga': 69.723,
		 'Ge': 72.61,
		 'As': 74.92159,
		 'Se': 78.96,
		 'Br': 79.904,
		 'Kr': 83.8,
		 'Rb': 85.4678,
		 'Sr': 87.62,
		 'Y': 88.90585,
		 'Zr': 91.224,
		 'Nb': 92.90638,
		 'Mo': 95.94,
		 'Tc': 97.9072,
		 'Ru': 101.07,
		 'Rh': 102.9055,
		 'Pd': 106.42,
		 'Ag': 107.8682,
		 'Cd': 112.411,
		 'In': 114.818,
		 'Sn': 118.71,
		 'Sb': 121.76,
		 'Te': 127.6,
		 'I': 126.90447,
		 'Xe': 131.29,
		 'Cs': 132.90543,
		 'Ba': 137.327,
		 'La': 138.9055,
		 'Ce': 140.115,
		 'Pr': 140.90765,
		 'Nd': 144.24,
		 'Pm': 144.9127,
		 'Sm': 150.36,
		 'Eu': 151.965,
		 'Gd': 157.25,
		 'Tb': 158.92534,
		 'Dy': 162.5,
		 'Ho': 164.93032,
		 'Er': 167.26,
		 'Tm': 168.93421,
		 'Yb': 173.04,
		 'Lu': 174.967,
		 'Hf': 178.49,
		 'Ta': 180.9479,
		 'W': 183.84,
		 'Re': 186.207,
		 'Os': 190.23,
		 'Ir': 192.22,
		 'Pt': 195.08,
		 'Au': 196.96654,
		 'Hg': 200.59,
		 'Tl': 204.3833,
		 'Pb': 207.2,
		 'Bi': 208.98037,
		 'Po': 208.9824,
		 'At': 209.9871,
		 'Rn': 222.0176,
		 'Fr': 223.0197,
		 'Ra': 226.0254,
		 'Ac': 227.0278,
		 'Th': 232.0381,
		 'Pa': 231.03588,
		 'U': 238.0289,
		 'Np': 237.048,
		 'Pu': 244.0642,
		 'Am': 243.0614,
		 'Cm': 247.0703,
		 'Bk': 247.0703,
		 'Cf': 251.0796,
		 'Es': 252.083,
		 'Fm': 257.0951,
		 'Md': 258.1,
		 'No': 259.1009,
		 'Lr': 262.11}
AMU.update({"foo": 0.0}) # Add this term for placeholder (TIP4P)