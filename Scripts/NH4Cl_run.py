""" NH4Cl Run file: Uses JBRQT potential and a constraint well""" 
import numpy as np
import initial_configuration as IC  # Generates molecules and initial state
import calculations as calc         # Calculates values from data
from sim import Simulation
import os

np.set_printoptions(suppress=True)

# Simulation save location (Don't need to write xyz, data_handler already adds it)
SAVE_FOLDER = "Sim" # Save folder in current directory
FILENAME = "NH4CL" # File save name
FILE_PATH = os.path.join(os.getcwd(), SAVE_FOLDER, FILENAME) # Save Location will be CurrentDirectory/saveFolder/fileName
SAVE = 10
CHECK = 1000

# Species List, currently using 2 NH4 and 2 Cl ions
speciesList =  {"NH4": 5,
                "Cl": 5}

# Create system and modify energy parameters 
Simulation.M = IC.Create_System.heterogeneous(speciesList, create_xyz = False)
class_list=  [calc.JBRQT, calc.N_Well]                  # Create potentials list to use
calc.N_Well.n = 10                                      # Change constraint n in (r/box)^2n 


sim_inputs = {
        "trials": 1000,                    # Number of trials
        "n_eq":0,                         # Equilibration trials
        "temperature":500,               # Temperature of run
        "save": SAVE,                     # No. trials before saving xyz files
        "step_size": 2.0,                 # Initial step size
        "step_check": 200,                # No. of trials before correcting step size
        "outputPath":FILE_PATH,           # Location of save file
        "seeding": True,                  # If True, will seed using given seed
        "debug": True,                    # If True, will print a log of moves and energies to help debug the code
        "potential_list": class_list,     # List of potentials (class names) to be used
        "diffCenter": False,              # Use different center point to calculate LJ 
        "LJCenter": 0                     # Atom to use in LJ distance calculations
        }

# Run system and store outputs 
Outputs = Simulation.run(**sim_inputs)   
list_of_outputs = Outputs.keys()
