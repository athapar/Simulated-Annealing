"""This program contains a class Simulation() which has the function run() to perform a Monte Carlo simualation

The run function has multiple inputs, of which trials is required and all else will take default values. 
The list of inputs and their descriptions are given in the sim_inputs dictionary

This dictionary can be copied to a run file and the inputs changed to suit the user's needs

It is advisable to keep the variable save high so that the act of saving data too often does 
not slow down the program. If the boolean "debug" is True the program will record all moves and
and create a log to help debug errors. This will slow down the program, however.  
"""

import numpy as np
import initial_configuration as IC  # Generates molecules and initial state
import parameters as pm             # Default parameters for simulation
from tqdm import tqdm
import data_handler as data         # Converts data into arrays
import calculations as calc         # Calculates values from data
from copy import deepcopy
import pandas as pd
import os
from copy import copy


#default simParameters for Homogeneous
N_PART = 2
DEFAULT_TRIALS = 1000
DEFAULT_T = 20 
DEFAULT_NEQ = 0 
DEFAULT_SAVE = 1000 
DEFAULT_STEP = 2.0 
DEFAULT_CHECK = 500
DEFAULT_SEED = 10 
DEFAULT_LJ_CENTER = 0
DEFAULT_POTENTIAL = [calc.Lennard_Jones]
DEFAULT_MOL_NAME = "H2O"
DEFAULT_FILENAME = "H2O" # File save name
DEFUALT_SAVE_FOLDER = "Sim"
DEFAULT_OUTPATH_PATH = os.path.join(os.getcwd(),DEFUALT_SAVE_FOLDER,DEFAULT_FILENAME)
DEFAULT_M_IN = []



class Simulation():
    
    
    
    homogeneous = True
    
    # Simulation Parameters 
    P_MT = 0.1               # Probability of translational mag-step
    P_RT = 0.1               # # Probability of rotational mag-step
    MAG_FACTOR = 10          # Factor by which to magnify translational step
    MAG_THETA = 180          # Maximum theta for rotational magstep range will be (-MAG_THETA, MAG_THETA)
    COL_DIST = 1.5           # Collision distance [A]
    N_PART = 2               # Number of particles by default (Used by IC.Create_System to make M)
    
    # Seeding random number generator (will skip if Anneal is True)
    USE_SEED = True                 # Use value for seed below if True, else use system time to seed
    SEED = DEFAULT_SEED             # Seed input (using system time and input needs to be integer)
    
    
    # Import parameters from parameteers file
    A = pm.TIP3P["A"]
    B = pm.TIP3P["B"]
    K_C= pm.K_C    # Coulomb constant
    SIGMA = pm.TIP3P["SIGMA"]
    EPSILON = pm.TIP3P["EPSILON"]
    K_B = pm.K_B # Boltzmann constant [kcal/mol-K]
    
        
    # Simulation parameters
    THETA = 90                # Max angle for non-magnified steps
    SAVE = 2                  # Number of iterations before saving
    CHECK = 100               # Check acceptance every IRATIO cycles
    BOX = 10
    
    # You can enter Simulation inputs for function directly into parenthesis or they can be unpacked from a dictionary using ** operator --> Simulation(**simInputs)
    # Showing dictionary here which is easier to copy and paste and edit in a different file
    # Values for N_part and trials are needed, while others take default keyword arguments (see inside parentheses of run())
    sim_inputs = {
            "trials": DEFAULT_TRIALS,                # Number of trials
            "temperature": DEFAULT_T,                # Temperature of run
            "n_eq": DEFAULT_NEQ,                     # Equilibration trials
            "save": DEFAULT_SAVE,                    # No. trials before saving xyz files
            "step_size": DEFAULT_STEP,               # Initial step size
            "step_check": DEFAULT_CHECK,             # No. of trials before correcting step size
            "outputPath": DEFAULT_OUTPATH_PATH,      # Location of save file
            "seeding": True,                         # If True, will perform random seed
            "debug": False,                          # If True, will print a log of moves and energies to help debug the code
            "potential_list": DEFAULT_POTENTIAL,     # List of potentials (class names) to be used
            "diffCenter": True,                      # Use different center point to calculate LJ 
            "LJCenter": DEFAULT_LJ_CENTER            # Atom to use in LJ distance calculations
            }
    
    M = DEFAULT_M_IN
        
    @classmethod 
    def update_M(cls,M_in):
        cls.M = deepcopy(M_in)              # This function updates M used by the simulation. Can also be done manually


    @classmethod
    def run(cls, trials, temperature= DEFAULT_T,  n_eq= DEFAULT_NEQ, save= DEFAULT_SAVE, step_size= DEFAULT_STEP, 
            step_check= DEFAULT_CHECK, outputPath= DEFAULT_OUTPATH_PATH, seeding= True, debug= False, 
            potential_list= DEFAULT_POTENTIAL, diffCenter= True,  LJCenter= DEFAULT_LJ_CENTER):           

        
        UP = calc.User_Potential  
        UP.class_list = potential_list  # Use potential with 3 terms
        potentials, potentials2 = UP.create_potentials()
        potential_names = UP.potential_names_list()
        
        # Seed using parameters above (seeding==False by default)
        if seeding:
            calc.seedRandom(cls.SEED, cls.USE_SEED)

        # Initialize arrays    
        Energy_Store=[]
        Accepted_Energies=[]
        dataLog = [] 
             
        # Import molecular array and calculate initial energy
        M = cls.M
        N_part = len(M)
        data.XYZ.create(M, path=outputPath)
        U_arr = UP.calculate_potentials(M, potentials)
        
        # Store the initial values
        storage_dict= data.create_storage_dict(U_arr,potential_names)
        Accepted_Energies.append(storage_dict)
        Energy_Store.append(storage_dict)
        

        # Initialize acceptance counters (format [accept, reject, total])
        total_acceptance_array = np.array([0,0,0])
        loop_acceptance_array = np.array([0,0,0])  
        count_iter = 0
        
        # k = trial iterator, i = molecule/atom iterator
        for k in tqdm(range(trials)):				
            for i in range(N_part):  
                count_iter += 1
                
                # Adjust step size
                if k%step_check == 0 and i == 0: 
                    step_size, acceptance_list = calc.update_step_size(step_size, temperature, loop_acceptance_array, step_check)
                
                
                # Find energies relative to molecule/atom i
                UI_arr = UP.calculate_potentials2(M, i, potentials2)
                
                # Save old position in case move is not accepted
                M[i].previous_position = deepcopy(M[i].position) # Use deepcopy to store previous position separately
               
    
                # Translation condition for magwalk
                if np.random.random() < cls.P_MT:
                    vector = cls.MAG_FACTOR*step_size * np.array([np.random.uniform(low=-1.0,high=1.0), np.random.uniform(low=-1.0,high=1.0), np.random.uniform(low=-1.0,high=1.0) ])
                    M[i].translate(vector)
                    M[i].update()
                else:
                    vector = step_size * np.array([np.random.uniform(low=-1.0,high=1.0), np.random.uniform(low=-1.0,high=1.0), np.random.uniform(low=-1.0,high=1.0)  ])
                    M[i].translate(vector)
                    M[i].update()
                    
                # Rotation condition (only applied to molecules, i.e. size>1)
                if M[i].size > 1:
                    if np.random.random() < cls.P_RT:
                        M[i].random_rotation(cls.MAG_THETA)
                    else:
                        M[i].random_rotation(cls.THETA)
                
                if k%save == 0:
                    data.XYZ.append(M, path = outputPath, homogeneous = cls.homogeneous)
    
    
                Collide = False
                
                # Loop to check if molecule collides with other molecules 
                for j in range(N_part):
                    if j == i:
                        pass
                    elif calc.distance(M[i].center, M[j].center) < 2*cls.COL_DIST:
                        Collide = True
                        break
                    else:
                        pass
    
                # REJECT MOVE BASED ON COLLISION
                if Collide:
    
                    if  k>= n_eq:
                        storage_dict= data.create_storage_dict(U_arr,potential_names)
                        Energy_Store.append(storage_dict)
                        
                        # Update counters
                        loop_acceptance_array += np.array([0,1,1])
                        total_acceptance_array += np.array([0,1,1])  

    
                    if debug:
                        log_dict = {"Trial":k, "Molecule Moved": i, "Random Number":"N/A", "Boltzmann":"N/A", "Result": "Reject", "Why":"Collision", "DeltaU":"N/A"}
                        dataLog.append({**log_dict, **storage_dict})
                  
                    M[i].position = M[i].previous_position
                    M[i].update()
    
                # If not colliding, perform energy test
                else:
                    # Find test energy and deltaU using the molecule/atom moved
                    test_UI_arr = UP.calculate_potentials2(M, i, potentials2) 
                    delta_UI_arr = test_UI_arr - UI_arr
                    test_U_arr = U_arr + delta_UI_arr 
    
                    rand_n = np.random.random()
                   # Set Boltzmann factor to 0 for T=0 so no high energy moves are accepted
                    if temperature == 0:
                        bolz = 0
                    else:
                        if isinstance(delta_UI_arr, np.ndarray):
                            bolz = np.exp(-delta_UI_arr[-1]/cls.K_B/temperature)
                        else: 
                            bolz = np.exp(-delta_UI_arr/cls.K_B/temperature)
                    
                    # Find deltaUI total depends on how long array of potential energies is
                    if isinstance(delta_UI_arr, np.ndarray):
                        delta_UI_tot = delta_UI_arr[-1]
                    else:
                        delta_UI_tot = delta_UI_arr
    
                
                    # Energy Decrease/No change = ACCEPT
                    if np.sign(delta_UI_tot) <= 0:
    
                        if k >= n_eq:                           
                            if k%save == 0:
                                data.XYZ.append(M, path=outputPath, homogeneous = cls.homogeneous)

                            # Update counters
                            loop_acceptance_array += np.array([1,0,1])
                            total_acceptance_array += np.array([1,0,1])    
                            
            
                            # Store new energies
                            storage_dict= data.create_storage_dict(test_U_arr,potential_names)
                            Accepted_Energies.append(storage_dict)
                            Energy_Store.append(storage_dict)
                            
                            
                            if debug:
                                log_dict = {"Trial":k, "Molecule Moved": i, "Random Number":"N/A", "Boltzmann":"N/A", "Result": "Accept", "Why":"Energy Decrease", "DeltaU":delta_UI_arr[-1]}
                                dataLog.append({**log_dict, **storage_dict})
    
                        # Update U when accepted
                        U_arr = copy(test_U_arr)
                        
    
                    # PASSES BOLTZMANN TEST = ACCEPT
                    elif rand_n <= bolz:

                        if k >= n_eq:
                            if k%save == 0:
                                data.XYZ.append(M, path=outputPath, homogeneous = cls.homogeneous)
                            
                            # Update counters
                            loop_acceptance_array += np.array([1,0,1])
                            total_acceptance_array += np.array([1,0,1])

                            # Store new energies
                            storage_dict= data.create_storage_dict(test_U_arr,potential_names)
                            Accepted_Energies.append(storage_dict)
                            Energy_Store.append(storage_dict)
                                            
                            if debug:
                                log_dict = {"Trial":k, "Molecule Moved": i, "Random Number":rand_n, "Boltzmann":bolz, "Result": "Accept", "Why":"Boltmann Test Pass", "DeltaU":delta_UI_arr[-1]}
                                dataLog.append({**log_dict, **storage_dict})
                                
                        # Update U when accepted
                        U_arr = copy(test_U_arr)
    
    
                    # Fails Boltzmann Test = REJECT
                    else:
    
                        if k >= n_eq:
                            if k%save == 0:
                                data.XYZ.append(M, path=outputPath, homogeneous = cls.homogeneous)
                            
                            # Update counters
                            loop_acceptance_array += np.array([0,1,1])
                            total_acceptance_array += np.array([0,1,1])

                            # Store new energies
                            rej_storage_dict= data.create_storage_dict(test_U_arr,potential_names)
                            Energy_Store.append(rej_storage_dict)
                           
                            
                            if debug:
                                storage_dict1 = data.create_storage_dict(U_arr,potential_names)
                                log_dict = {"Trial":k, "Molecule Moved": i, "Random Number":rand_n, "Boltzmann":bolz, "Result": "Reject", "Why":"Boltmann Test Fail", "DeltaU":delta_UI_arr[-1]}
                                dataLog.append({**log_dict, **storage_dict1})
    
  
                        M[i].position = M[i].previous_position
                        M[i].update()
        
        # Allow variable no. of outputs using dictionary, have more outputs when debugging
        if debug: 
            Outputs = {"Energy_Store":pd.DataFrame(Energy_Store), 
                       "Accepted_Energies": pd.DataFrame(Accepted_Energies),
                       "Acceptance": total_acceptance_array,
                       "M_fin": M,
                       "step_size": step_size,
                       "dataLog": pd.DataFrame(dataLog[1:])}
        else:
            Outputs = {"Energy_Store":pd.DataFrame(Energy_Store),
                       "Accepted_Energies": pd.DataFrame(Accepted_Energies),
                       "Acceptance": total_acceptance_array,
                       "M_fin": M,
                       "step_size": step_size}
             
        return Outputs  


# Simulation.M = IC.Create_System.homogeneous(N_PART, moleculeName="H2O", Mol_Class=IC.Molecule, create_xyz=False)
# Outputs = Simulation.run(**Simulation.sim_inputs)

