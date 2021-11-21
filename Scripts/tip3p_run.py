""" TIP3P Run file with special class added for model 3-site model""" 

import numpy as np
import initial_configuration as IC  # Generates molecules and initial state
import parameters as pm             # Default parameters for simulation
import calculations as calc         # Calculates values from data
from copy import deepcopy
from sim import Simulation
import os


# Sim save location (Don't need to write xyz, function already adds it)
MOL_NAME = "H2O" # This is fed to molecule class or custom class to create instances
SAVE_FOLDER = "Sim" # Save folder in current directory
FILENAME = "TIP3P" # File save name
FILE_PATH = os.path.join(os.getcwd(), SAVE_FOLDER, FILENAME) # Save Location will be CurrentDirectory/saveFolder/fileName
SAVE = 1
CHECK = 1000
N_PART = 2

# Import TIP3P paramters from paramters file
SIGMA = pm.TIP3P["SIGMA"]
EPSILON = pm.TIP3P["EPSILON"]



class_list = [calc.Lennard_Jones, calc.Coulombic, calc.N_Well]  # Use 3 potentials
calc.Lennard_Jones.sigma = SIGMA                 # Modify sigma and epsilon to match TIP3P
calc.Lennard_Jones.epsilon = EPSILON 
calc.Lennard_Jones.diff_center = True           # Don't use COM to calculate r, 
calc.Lennard_Jones.center_atom = 0              # Use oxygen to calculate distance for LJ 6-12
calc.N_Well.n = 2                               # Use quadwell potential (r/C)^2n = (r/C)^4

### MODIFIED MOLECULE CLASS FOR SIMULATION
# Create modified molecule class
class Water_Mol(IC.Molecule):
    # Modify __init__() to include charges and specify oxygen and hydrogen positions separately
    def __init__(self, name = "H2O"): 
        super().__init__("H2O")
        self.name = name

        # Add partial charges based on TIP3P model
        self.charges = [pm.TIP3P["Q_O"], pm.TIP3P["Q_H"], pm.TIP3P["Q_H"]] # Arrangement is [O H H]
        
        # Test position for particle
        self.previous_position = deepcopy(self.position)

    def update(self):
        super().update()

    def random_rotation(self, angle):
        x_angle = angle*np.random.uniform(low=-1.0,high=1.0)
        y_angle = angle*np.random.uniform(low=-1.0,high=1.0)
        z_angle = angle*np.random.uniform(low=-1.0,high=1.0)


        self.xrotation(x_angle)
        self.update()

        self.yrotation(y_angle)
        self.update()

        self.zrotation(z_angle)
        self.update()

        return self.position
### END MODIFIED CLASS


# Create M
Simulation.M = IC.Create_System.homogeneous(N_PART, moleculeName="H2O", Mol_Class=Water_Mol, create_xyz=False)


# You can enter Simulation inputs for function directly into parenthesis or they can be unpacked from a dictionary using ** operator --> Simulation(**simInputs)
# Showing dictionary here to help to explain the individual inputs more easily
simInputs = {
        "trials": 100,                    # Number of trials
        "n_eq":0,                         # Equilibration trials
        "temperature":100,                # Temperature of run
        "save": SAVE,                     # No. trials before saving xyz files
        "step_size": 0.5,                 # Initial step size
        "step_check": 100,                # No. of trials before correcting step size
        "outputPath":FILE_PATH,           # Location of save file
        "seeding": True,                  # If True, will seed using given number
        "debug": False,                   # If True, will print a log of moves and energies to help debug the code
        "potential_list": class_list,     # List of potential classes to be used
        "diffCenter": True,               # Use different center point to calculate LJ 
        "LJCenter": 0                     # Atom to use in LJ distance calculations
        }

# Run and store outputs to a dictionary variable
Outputs = Simulation.run(**simInputs)    
M = Outputs["M_fin"]


list_of_outputs = Outputs.keys()