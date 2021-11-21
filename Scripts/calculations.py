"""
This Program contains functions used for calculations in both the setup and during
the course of a simulation. Most of the functions in this library require the input of
a molecular system, M, to iterate on.
"""
import numpy as np
import parameters as pm
import time


# Takes any two position vectors and finds distance between them
def distance(v1,v2):
    diffv = v1-v2
    D = 0
    for iter1,term in enumerate(diffv):
        D += term**2
    D = np.sqrt(D)
    return D


# Finds distance between molecules using COM
def object_distance(M1,M2):
    return distance(M1.center, M2.center)


def seedRandom(inputSeed,condition):
    if condition:
        np.random.seed(inputSeed)
    else:
        np.random.seed(int(time.time()))
    return None


def system_center(M):
    array_sum = np.array([0.0, 0.0, 0.0])
    total_mass = 0.0 
    
    for obj in M:
        array_sum += obj.mass * obj.center 
        total_mass += obj.mass
        
    return array_sum/total_mass
    
def update_step_size(t_step, temp, acceptance_list, step_check_num):
    """ Takes inputs step-size, temperature, acceptance_list, step_check"""
    accepted,rejected,total = acceptance_list
    if total == 0:
        pass
    else:
        acceptance = accepted/total

        if acceptance <= 0.5:
            if temp == 0:
                t_step = 0.01
            elif acceptance <= 0.01:
                if t_step < 0.001:
                    pass
                else:
                    t_step *= 0.1
            else:
                t_step *= 0.8 
                
            acceptance_list = np.array([0,0,0])

        elif acceptance > 0.5:
            if acceptance >= 0.7:
                t_step *= 5
            else:
                t_step *= 1.2                
            acceptance_list = ([0,0,0])
            
    return t_step, acceptance_list



############################################################################
# Potential Classes
############################################################################


class Lennard_Jones():
    """ Calculates Lennard Jones 6-12 potential of the form U = 4*epsilon((sigma/r)^12 - (sigma/r)^6) for a system of molecules"""
    
    # Class variables used as parameters by class methods using syntax cls.sigma, etc.
    sigma = pm.SIGMA                 # Parameter sigma for LJ 6-12
    epsilon = pm.EPSILON             # Parameter epsilon for LJ 6-12
    calculate_atom_by_atom = False   # If True, calculate potential using all the atoms, not molecules only
    diff_center = False              # If molecule only: choose if distance r calculated using default COM of molecule or particular atom
    center_atom = 0                  # if particular atom, give atom position in atom_names list
    print_string = "U_LJ"            # Used when printing outputs

    #########################
    # Potential 1 (Entire system)
    #########################    
    
    # class method decorator allows function to use class variables ^^^
    @classmethod 
    def potential(cls, M): # Will be in kcal/mol
        En = 0
    
        if cls.calculate_atom_by_atom == False:
            for iter1, obj1 in enumerate(M[:-1]):
                for iter2,obj2 in enumerate(M[iter1+1:], start = iter1+1):
                    if cls.diff_center:
                        C1 = obj1.position[cls.center_atom]
                        C2 = obj2.position[cls.center_atom]
                    else:
                        C1 = obj1.center
                        C2 = obj2.center
                        
                    r =  distance(C1,C2) # Currently takes first molecule as center
                    r6, s6 = r**6, cls.sigma**6
                    r12, s12 = r6**2, s6**2
                    En += s12/r12 - s6/r6
                    
            En *= 4 * cls.epsilon
           
            

        else:
            for iter1, obj1 in enumerate(M[:-1]):
                for iter2,obj2 in enumerate(M[iter1+1:], start = iter1+1):
                    for a in range(obj1.size):
                        for b in range(obj2.size):
                            
                            if cls.epsilon[a][b]==0 or cls.sigma[a][b]==0:
                                pass
                            
                            else:
                                dist = distance(obj1.position[a],obj2.position[b])
                                r6, s6 = dist**6, cls.sigma[a][b]**6
                                r12, s12 = r6**2, s6**2
                                En += 4*cls.epsilon[a][b]*(s12/r12 - s6/r6)
            
        return En



    #########################
    # Potential 2 (Relative to one molecule only --> fast computations)
    #########################  
    
    # Used to compute deltaU quickly if object i was moved (computes U relative to i only)
    @classmethod # Class method decorator
    def potential2(cls, M, i): # Will be in kcal/mol
        En = 0
        if cls.calculate_atom_by_atom == False:
            
            for iter1, obj in enumerate(M): 
                if i == iter1:
                    pass
                
                else:
                    if cls.diff_center: # Pick atom to calculate distances
                        C1 = M[i].position[cls.center_atom]
                        C2 = obj.position[cls.center_atom]
                    else:
                        C1 = M[i].center
                        C2 = obj.center
                    
                    dist = distance(C1,C2) # Currently takes first molecule as center
                    r6, s6 = dist**6, cls.sigma**6
                    r12, s12 = r6**2, s6**2
                    En += s12/r12  - s6/r6
                    
            En *= 4 * cls.epsilon
                    
    
        else:
            
            for iter1,obj in enumerate(M):
                if i == iter1: 
                    pass
            
                else:
                    for a in range(M[i].size):
                        for b in range(obj.size):
                            if cls.epsilon[a][b] == 0 or cls.sigma[a][b] == 0:
                                pass
                            else:
                                dist = distance(M[i].position[a], obj.position[b])
                                r6, s6 = dist**6, cls.sigma[a][b]**6
                                r12, s12 = r6**2, s6**2
                                En += 4*cls.epsilon[a][b]*(s12/r12 - s6/r6)
                                
    
        return En


class Coulombic():
    """ Coulombic potential: calculates potential in the form U = k_c*qi*qj/r_ij for all atoms in a system"""
    
    # Class variables used as parameters by class methods using cls.k_c, etc.
    k_c = pm.K_C                      # Coulomb constant (kcal/mol)
    print_string = "U_C"              # Used when printing outputs
    
    #########################
    # Potential 1 (Entire System)
    #########################       
    
    # class method decorator allows function to use class variables ^^^
    @classmethod
    def potential(cls,M):
        En = 0
        for iter1, obj1 in enumerate(M[:-1]):
            for iter2,obj2 in enumerate(M[iter1+1:],start=iter1+1):
                EIJ = 0
                for a in range(obj1.size):
                    for b in range(obj2.size):
                        dist = distance(obj1.position[a],obj2.position[b])
                        
                        EIJ += obj1.charges[a] * obj2.charges[b] / dist
                
                EIJ *= cls.k_c
                En += EIJ
        
        return En

    #########################
    # Potential 2 (Relative to one molecule only --> fast computations)
    #########################   

    @classmethod
    def potential2(cls, M, i):
        En = 0
        for iter1,obj in enumerate(M):
            EIJ = 0
            if iter1 == i:
                pass
            else:
                for a in range(M[i].size):
                    for b in range(obj.size):
                        dist = distance(M[i].position[a], obj.position[b])
                        EIJ += M[i].charges[a] * obj.charges[b] / dist
    
                EIJ *= cls.k_c # Coulombic Term
                En += EIJ
        return En
    
    
class N_Well():
    """ Constraint potential in the form U = r^(2n) to a reference point"""
    # Function parameters
    n = 2                             # input exponent for potential
    box = 10                          # Constraint box size
    ref = np.zeros(3)                 # Reference point for distance (can be origin, system COM, etc.)
    print_string = "U_Constr"         # Used when printing outputs 
    
    #########################
    # Potential 1
    #########################   
    @classmethod
    def potential(cls, M):
        En = 0 
        for obj in M:
            En += (distance(cls.ref, obj.center) / cls.box)**(2*cls.n)
        return En

    #########################
    # Potential 2
    #########################  
    @classmethod 
    def potential2(cls, M, i):
        En = (distance(cls.ref, M[i].center) / cls.box)**(2*cls.n)
        return En
    
    
    #########################
    # Description of inputs (not used in calculations anymore)
    ######################### 
    
    inputs = {"M":[],                   # Leave blank for input, program will automatically update
              "box": 10}                # Box size for simulation
    
    inputs2 = {"M":[],                    # Leave blank for input, program will automatically update
               "i": 0,                    # Molecule number (loop will update this on its own)
               "box": 10}                 # Box size for simulation
  
    
class TIP3P():
    """ Calculate potential via equation: U = k_c*qi*qj/r + A/rOO^12 - B/rOO^6"""
    k_c = pm.K_C             # Coulomb constant
    A = pm.TIP3P["A"]        # Parameter A for TIP3P 
    B = pm.TIP3P["B"]        # Parameter B for TIP3P potential
    ref = np.zeros(3)
    
    @classmethod
    def compute(cls,M): # Will be in kcal/mol
        E = 0
        for i in range(len(M)-1):
            for j in range(i+1,len(M)):
    
                # Inter Oxygen distance
                r_oo = distance(M[i].Ox, M[j].Ox)
                r6 = r_oo ** 6
                r12 = r6 ** 2
    
    
                # Initialize variable for charge potential energy
                EIJ = 0
                for a in range(3):
                    for b in range(3):
                        EIJ += M[i].charges[a] * M[j].charges[b] / distance(M[i].position[a],M[j].position[b])
    
                EIJ = cls.k_c * EIJ # Coulombic Term
                EIJ = EIJ + cls.A/r12 - cls.B/r6	 # Pair Potential Energy between
                E += EIJ 
            E += distance(np.zeros(3),M[i].center) # Add Pair Potential to Total Potential Energy
        E += distance(np.zeros(3),M[-1].center)
        return E    
    
    @classmethod
    def compute2(cls,M,i):
        En = 0
        for j in range(len(M)):
            if j == i:
                pass
            else:
                EIJ = 0
    
                r_oo = distance(M[i].center, M[j].center)
                r6, s6 = r_oo**6, cls.sigma**6
                r12, s12 = r6**2, s6**2
                r6 = r_oo ** 6
                r12 = r6 ** 2
                for a in range(3):
                    for b in range(3):
                        EIJ += M[i].charges[a] * M[j].charges[b] / distance(M[i].position[a],M[j].position[b])
                EIJ = cls.k_c * EIJ
                EIJ += 4*cls.epsilon*(s12/r12 - s6/r6)
                En += EIJ			
        En += distance(np.zeros(3),M[i].center)
        return En
    
    
class User_Potential():
    """ Creates custom potential based on user input: By default calculates using LJ, Coulombic, N_Well"""
    class_list  = [Lennard_Jones, Coulombic, N_Well]
    sum_at_end = True
          
    @classmethod
    def create_potentials(cls):
        f1 = []
        f2 = []
        
        for potential_class in cls.class_list:
            f1.append(potential_class.potential)
            f2.append(potential_class.potential2)
        return f1, f2 
    
    @classmethod
    def calculate_potentials(cls, M, funcs):
        if len(funcs) > 1: 
            U_list = []
            for func in funcs:
                result = func(M)
                if type(result) == list or type(result) == np.ndarray:
                    U_list = [*U_list, *result]
                else:
                    U_list.append(func(M))
            if cls.sum_at_end == True:
                U_list.append(sum(U_list))     
        else:
            U_list = funcs[0](M)
        return np.array(U_list)
    
    @classmethod      
    def calculate_potentials2(cls, M, i, funcs):
        if len(funcs) > 1:
            U_list = []
            for func in funcs:
                result = func(M, i)
                if type(result) == list or type(result) == np.ndarray:
                    U_list = [*U_list, *result]
                else:
                    U_list.append(func(M,i))
            if cls.sum_at_end == True:
                U_list.append(sum(U_list))
        else:
            U_list = funcs[0](M,i)
        return np.array(U_list)
    
    @classmethod
    def potential_names_list(cls):
        if len(cls.class_list) > 1:
            names_list = []
            for pot in cls.class_list:
                if type(pot.print_string) == list:
                    names_list = [*names_list, *pot.print_string]
                else:
                    names_list.append(pot.print_string)
                    
            if cls.sum_at_end == True:
                names_list.append("U_tot")
        else:
            names_list = cls.class_list[0].print_string
        return names_list 


class JBRQT():
    A, B, C, D, Q = pm.AmmCl["A"], pm.AmmCl["B"], pm.AmmCl["C"], pm.AmmCl["D"], pm.AmmCl["Q"]
    k_c = pm.K_C
    print_string = ["U_buck", "U_coul", "ULJ6", "ULJ12", "U_tot"]
    test_val = 10
    

    #########################
    # Potential 1
    #########################     
    
    @classmethod
    def potential(cls, M): 
        U_buck = 0.0
        U_coul = 0.0 
        ULJ6 = 0.0
        ULJ12 = 0.0      
        for iter1, mol1 in enumerate(M[:-1]): # Using iter1 and iter2 so that i in loop is not changed
            for iter2, mol2 in enumerate(M[iter1+1:],start=iter1+1):        
                for a in range(mol1.size):
                    for b in range(mol2.size):
                        if mol1.size == 1:
                            if mol2.size == 1:
                                dist = distance(mol1.position, mol2.position)
                            else:
                                dist = distance(mol1.position, mol2.position[b])
                        elif mol2.size == 1:
                            dist = distance(mol1.position[a], mol2.position)
                        else:
                            dist = distance(mol1.position[a], mol2.position[b])
                            
                        r6 = dist**6
                        r12 = r6**2
                        k1 = mol1.atom_names[a]  # atom1 name or key
                        k2 = mol2.atom_names[b]  # atom2 name or key
                        U_buck += cls.A[k1][k2] * np.exp(-cls.B[k1][k2]*dist)
                        U_coul += cls.k_c*cls.Q[k1][k2] / dist
                        ULJ6 += -cls.C[k1][k2]/r6 
                        ULJ12 += cls.D[k1][k2]/r12
                        
                
                
        U_arr = [U_buck, U_coul, ULJ6, ULJ12] # Put all energies into array
        U_arr.append(sum(U_arr))    # Add total energy to array
        return np.array(U_arr)



    @classmethod
    def potential2(cls, M, i): 
        U_buck = 0.0
        U_coul = 0.0 
        ULJ6 = 0.0
        ULJ12 = 0.0   
        for iter1, obj in enumerate(M): # Using iter1 and iter2 so that i in loop is not changed
            if iter1 == i:
                pass
            else:     
                for a in range(M[i].size):
                    for b in range(obj.size):
                        if M[i].size == 1:
                            if obj.size == 1:
                                dist = distance(M[i].position, obj.position)
                            else:
                                dist = distance(M[i].position, obj.position[b])
                        elif obj.size == 1:
                            dist = distance(M[i].position[a], obj.position)
                        else:
                            dist = distance(M[i].position[a], obj.position[b])
                        
                        
                        
                        r6 = dist**6
                        r12 = r6**2
                        k1 = M[i].atom_names[a]  # atom1 name or key
                        k2 = obj.atom_names[b]  # atom2 name or key
                        
                        U_buck += cls.A[k1][k2] * np.exp(-cls.B[k1][k2]* dist)
                        U_coul += cls.k_c*cls.Q[k1][k2] / dist
                        ULJ6 += -cls.C[k1][k2]/r6 
                        ULJ12 += cls.D[k1][k2]/r12
                        
                        
        U_arr = [U_buck, U_coul, ULJ6, ULJ12] # Put all energies into array
        U_arr.append(sum(U_arr))    # Add total energy to array
        return np.array(U_arr)

