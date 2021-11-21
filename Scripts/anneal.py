""" This file contains the annealing program: the default data is currently for a homogeneous simulation
of TIP3P water. The potentials are also adjusted for TIP3P water. This file can be modified for another homogeneous
simulation or made into a heterogenous simulation, though the NH4Cl_anneal.py may be a better starting point. """

import data_handler as data
import numpy as np
import os
import calculations as calc
import pandas as pd 
import tip3p_run as t3
import initial_configuration as IC
import parameters as pm
from sim import Simulation

calc.seedRandom(10,True)
# Folder names
SAVE_FOLDER = "TIP3P" # <---Folder name that may need to be changed (Currently saved inside Sim)
MAIN_FOLDER = "Sim" # Usually the folder where other sims are already saved
DATA_FOLDER = "Data" 
XYZ_FOLDER = "Movie"
PLOT_FOLDER = "Plots"
# Create file path names (CWD = current working directory)
TRIAL_PATH = os.path.join(os.getcwd(), MAIN_FOLDER, SAVE_FOLDER) # CWD/mainFolder/saveFolder
DATA_PATH = os.path.join(TRIAL_PATH, DATA_FOLDER) # CWD/mainFolder/saveFolder/dataFolder
XYZ_PATH = os.path.join(TRIAL_PATH, XYZ_FOLDER) # CWD/mainFolder/saveFolder/XYZFolder
PLOT_PATH = os.path.join(TRIAL_PATH, PLOT_FOLDER) # CWD/mainFolder/saveFolder/plotFolder

paths = [TRIAL_PATH, DATA_PATH, XYZ_PATH, PLOT_PATH]
data.createPaths(paths) # Creates paths, or ignores if they already exist


RUNS = 100000
NEQ = 10000
SAVE = 1000
CHECK = 300
DEFAULT_STEP = 2.0
showFigs = False 
Tmax= 1000 # K
N_PT = 5# Temperatures per tooth
N_T = 2# Number of teeth 
alpha  = 0.7 # Factor by which max T decreases in subsequent teeth
N_PART = 2 
SAVE= 100



Water_Mol= t3.Water_Mol
class_list = [calc.Lennard_Jones, calc.Coulombic, calc.N_Well]
calc.Lennard_Jones.sigma = pm.TIP3P["SIGMA"]
calc.Lennard_Jones.epsilon = pm.TIP3P["EPSILON"]
calc.N_Well.n = 10
Simulation.M = IC.Create_System.homogeneous(N_PART, moleculeName="H2O", Mol_Class=Water_Mol, create_xyz=False)



class Anneal():
   
    
    sim_inputs = {
        "trials": RUNS,                   # Number of trials
        "n_eq": NEQ,                      # Equilibration trials
        "temperature":1000,               # Temperature of run
        "save": SAVE,                     # No. trials before saving xyz files
        "step_size": DEFAULT_STEP,        # Initial step size
        "step_check": CHECK,              # No. of trials before correcting step size
        "outputPath":TRIAL_PATH,          # Location of save file
        "seeding": True,                  # If True, will perform random seed
        "debug": False,                   # If True, will print a log of moves and energies to help debug the code
        "potential_list": class_list,     # List of potential (classes) to be used
        "diffCenter": True,               # Use different center point to calculate LJ 
        "LJCenter": 0                     # Atom to use in LJ distance calculations
        }
    
    
    @classmethod
    def run(cls):
        
        dataStore = []
        temperatures = []
        T1 = Tmax
        
        for P in range(1, N_T+1):
            T = T1
            delT = T1/(N_PT - 1)
        
            for N in range(N_PT):
                print("\n\nP=",P,"N=",N, "T=",T)
                PNTString = "P"+ str(P) + "_N" + str(N) + "_T" + str(T)  # Use in file names to identify N, P, Temperature
                #Set up xyz filename
                xyz_file = os.path.join(XYZ_PATH,PNTString)   
        
                if P == 1 and N == 0:
                    cls.sim_inputs["temperature"] = T
                    cls.sim_inputs["outputPath"] = xyz_file
                    # temperatures.append(T)
                    
                    
                    # RUN SIMULATION
                    Outputs = Simulation.run(**cls.sim_inputs)
                    M_fin = Outputs["M_fin"]
                    # dataStore.append([T, M_fin])
                    
                    U_accepted = Outputs["Accepted_Energies"].iloc[:,-1]
                    all_accepted_U = U_accepted # Save to a running tally
                    print("\n", U_accepted)
                    print(Outputs["Accepted_Energies"])
                    

                    
                    # Create Movie file for storing final position
                    data.XYZ.create(M_fin,path=os.path.join(XYZ_PATH,"Movie"))
                
                else:
                    if np.round(T) == 0:
                        T = 0
                        
                    # Simulation.M = deepcopy(M_fin )
                    cls.sim_inputs["temperature"] = T
                    cls.sim_inputs["outputPath"] = xyz_file
                    # temperatures.append(T)
                    
                    
                    # RUN SIMULATION
                    Outputs = Simulation.run(**cls.sim_inputs)
                    M_fin = Outputs["M_fin"]
                    # dataStore.append([T, M_fin])
                    
                    U_accepted = Outputs["Accepted_Energies"].iloc[:,-1]
                    all_accepted_U = pd.concat([all_accepted_U , U_accepted])
                    
                    # Create Movie file for storing final position
                    data.XYZ.append(M_fin,path=os.path.join(XYZ_PATH,"Movie"))
                  
                    
                # Save all and accepted energies: LJ, Coulombic, Constraint, Total
                Outputs["Energy_Store"].to_csv(os.path.join(DATA_PATH,"Energy_Store_"+PNTString+".csv"))
                Outputs["Accepted_Energies"].to_csv(os.path.join(DATA_PATH,"Accepted_Energies_"+PNTString+".csv"))
                temperatures.append(T)
                dataStore.append([T, M_fin])
                    
                if T == 0:
                    data.XYZ.create(M_fin, path= os.path.join(XYZ_PATH,"Final_Structure_P={}".format(P)))
                else:
                    T -= delT
                
            T1 = alpha*T1
            
        output_dict = {"temperatures": temperatures, "data_store": dataStore, "M_fin": M_fin, "all_accepted_U": all_accepted_U }
        
        return output_dict
                
             
        
results = Anneal.run()


