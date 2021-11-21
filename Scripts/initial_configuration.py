""" 
Initial Array Generator

Creates molecule and atom objects and stores them into an array (M) 
Each object has associated properties including mass, atom_names list, and position array

Position Array Shape (object composed of m atoms):
    x   y   z
1:[x1  y1  z1]
2:[x2  y2  z2]
...
j:[xj yj  zj]
...
m:[xm  ym  zm]

Particle positions in array:
M[i] = object i 
M[i].position = positions for all atoms in object i
M[i].position[j] = will select row j in position array: [xj, yj, zj] for atom j in object i

"""

## Import modules and define parameters
import numpy as np
import parameters as pm
import math
import calculations as calc
import data_handler as data
from copy import deepcopy
import os


MAXT = 0.2 # Maximum translation
MAXR = 360.0 # Maximum rotation
COL_DIST = 2 # [A] Collision distance
MAX_SEARCH = 100 # Max number of attempts to find a position that isn't colliding with other particles
INCREASE_MOVE = 1.1 # Factor to increase move size if algorithm fails to find position without colliding

# Parameters for creating save file fed to create_system function (Don't need to write xyz, function already adds it)
molName = "H2O" # This is fed to molecule class or custom class to create instance
saveFolder = "Sim" # Save folder in current directory
fileName = "mixed_sim" # File save name
filePath = os.path.join(os.getcwd(),saveFolder,fileName) # Save Location will be CurrentDirectory/saveFolder/fileName
customPath = os.path.join(os.getcwd(),saveFolder,"Custom.xyz") # Path of custom file



class Atom():
    # Initial Parameters
    def __init__(self, name, custom_path="Sim/Custom.xyz"):
        self.name = name # Specify type of molecule by name
        if self.name.upper() == "CL" or self.name.upper() == "CHLORINE":
            # Basic parameters
            self.name = "Cl"
            self.size = 1 						           # Number of atoms
            self.atom_names = ["Cl"] 	                   # Keys for atoms
            self.col_sphere = pm.Collision_Sphere["Cl"]    # Collision-sphere [A]
            self.mass = pm.AMU["Cl"]                       # Molecular mass
            self.center = [0.0,0.0,0.0]                    # Initial center is origin

            self.position = np.array([0.0, 0.0, 0.0])      # current in the form [A1 A2 A3...] where Ai = [xi yi zi]^T


            self.previous_position = self.position.copy()
            self.type = "atom"
            
        elif self.name.upper() == "F" or self.name.upper() == "FLUORINE":
            # Basic parameters
            self.name = "F"
            self.size = 1 						           # Number of atoms
            self.atom_names = ["F"] 	                   # Keys for atoms
            self.col_sphere = pm.Collision_Sphere["F"] 	   # Collision-sphere [A]
            self.mass = pm.AMU["F"]                        # Molecular mass
            self.center = [0.0,0.0,0.0]                    # Initial center is origin

            self.position = np.array([0.0, 0.0, 0.0])      # current in the form [A1 A2 A3...] where Ai = [xi yi zi]^T


            self.previous_position = self.position.copy()
            self.type = "atom"

        elif self.name.upper() == "BR" or self.name.upper() == "BROMINE":
            # Basic parameters
            self.name = "Br"
            self.size = 1                                  # Number of atoms
            self.atom_names = ["Br"]                       # Keys for atoms
            self.col_sphere = pm.Collision_Sphere["Br"]    # Collision-sphere [A]
            self.mass = pm.AMU["Br"]                       # Molecular mass
            self.center = [0.0,0.0,0.0]                    # Initial center is origin

            self.position = np.array([0.0, 0.0, 0.0])      # current in the form [A1 A2 A3...] where Ai = [xi yi zi]^T


            self.previous_position = self.position.copy()
            self.type = "atom"


        elif self.name.upper() == "I" or self.name.upper() == "IODINE":
            # Basic parameters
            self.name = "I"
            self.size = 1                                  # Number of atoms
            self.atom_names = ["I"]                       # Keys for atoms
            self.col_sphere = pm.Collision_Sphere["I"]    # Collision-sphere [A]
            self.mass = pm.AMU["I"]                       # Molecular mass
            self.center = [0.0,0.0,0.0]                    # Initial center is origin

            self.position = np.array([0.0, 0.0, 0.0])      # current in the form [A1 A2 A3...] where Ai = [xi yi zi]^T


            self.previous_position = self.position.copy()
            self.type = "atom"

        else:
            self.type = None

    def translate(self,delta):
        self.position += np.array(delta)

    def update(self):
        self.center = self.position


    # Copies previous position and then randomly moves the particle
    def random_move(self, dist=1, theta=90):
        self.translate([dist*(2*np.random.random()-1),  # Translate each randomly
                        dist*(2*np.random.random()-1),
                        dist*(2*np.random.random()-1)]) 		
        self.update()


# Creates a molecule with initial position given by self.position under the __init__() function
# Each molecule has functions to find COM and positions with COM as origin(COM_FRAME()) for rotations
class Molecule():
    # Initial Parameters
    def __init__(self, name, custom_path=customPath):
        self.name = name # Specify type of molecule by name
        if self.name.upper() == "H2O" or self.name.upper() == "WATER":
            # Basic parameters
            self.name = "H2O"
            self.size = 3 						# Number of atoms
            self.atom_names = ["O", "H", "H"] 	    # Keys for atoms
            self.col_sphere = pm.Collision_Sphere["H2O"] 	# Collision-sphere [A]
            self.mass = 2*pm.AMU["H"] + pm.AMU["O"] # Molecular mass
            self.center = [0,0,0] # Initial center is origin

            self.position = np.array([[0,0.06556811, 0], #  Default position
                                      [-0.757,-0.52043189, 0],
                                      [0.757, -0.52043189, 0]])

            # The molecule is initially centered at zero so the com_frame is the position matrix
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        elif self.name.upper() == "SO2":
            self.name = "SO2"
            # Basic Parameters
            self.size = 3  			# Number of atoms
            self.atom_names = ["S", "O", "O"]	# Keys for atoms
            self.col_sphere = pm.Collision_Sphere["SO2"]     # Collision-sphere [A]
            self.mass = pm.AMU["S"]+  2*pm.AMU["O"]
            self.position = np.array([[0,0,0],
                                      [0,1.2371,0.7215],
                                      [0,-1.2371,0.7215]]) # Initialize position
            
            self.center = [0,0,0.3603716]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()			
            self.type = "molecule"


        elif self.name.upper() == "N2" or self.name.upper() == "NITROGEN":
            self.name = "N2"
            # Basic Parameters
            self.size = 2  			# Number of atoms
            self.atom_names = ["N", "N"]	# Keys for atoms
            self.col_sphere = pm.Collision_Sphere["N2"]       # Collision-sphere [A]
            self.mass = 2*pm.AMU["N"]
            self.position = np.array([[0,0.54875,0],
                                      [0,-0.54875,0]]) # Initialize position
            self.center = [0,0,0]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"

    
        elif self.name.upper() == "HCL":
            self.name = "HCl"
            # Basic Parameters
            self.size = 2  			# Number of atoms
            self.atom_names = ["Cl", "H"]	# Keys for atoms
            self.col_sphere = pm.Collision_Sphere["HCl"]       # Collision-sphere [A]
            self.mass = 2*pm.AMU["N"]
            self.position = np.array([[0, 0, 0],
                                      [0, 0, 1.2746]]) # Initialize position
            self.center = [0,0,0]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()		
            self.type = "molecule"


        elif self.name.upper() == "NH3" or self.name.upper() == "AMMONIA":
            self.name = "NH3"

            # Basic Parameters
            self.size = 4
            self.atom_names = ["N", "H", "H", "H"]
            self.col_sphere = pm.Collision_Sphere["NH3"]
            self.mass = pm.AMU["N"] + 3*pm.AMU["H"]

            self.position = np.array([[0, 0, 0.1111],
                                      [0, 0.9316, -0.2592],
                                      [0.8068, -0.4658, -0.2592],
                                      [-0.8068, -0.4658, -0.2592]])

            self.center = [0,0,0.04535749]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        elif self.name.upper() == "NH4" or self.name.upper() == "AMMONIUM":
            self.name = "NH4"

            # Basic Parameters
            self.size = 5
            self.atom_names = ["N", "H", "H", "H", "H"]
            self.col_sphere = pm.Collision_Sphere["NH4"]
            self.mass = pm.AMU["N"] + 4*pm.AMU["H"]

            self.position = np.array([[0,0,0],
                                      [0.5939,0.5939,0.5939],
                                      [-0.5939,-0.5939,0.5939],
                                      [-0.5939,0.5939,-0.5939],
                                      [0.5939,-0.5939,-0.5939]])

            self.center = [0,0,0]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        elif self.name.upper() == "SF4":
            self.size = 5
            # Basic Parameters
            self.atom_names = ["S", "F", "F", "F", "F"]
            self.col_sphere = pm.Collision_Sphere["SF4"]
            self.mass = pm.AMU["S"] + 4*pm.AMU["F"]
            self.position = np.array([ [0,0,0.3825],
                                      [0,1.6255,0.2401],
                                      [0,-1.6255,0.2401],
                                      [1.2055,0,-0.581],
                                      [-1.2055,0,-0.581] ])
            self.center = [0,0,-0.00606887]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        elif self.name.upper() == "SF6":
            self.size = 7
            # Basic Parameters
            self.atom_names = ["S", "F", "F", "F", "F", "F", "F"]
            self.col_sphere = pm.Collision_Sphere["SF6"]
            self.mass = pm.AMU["S"] + 6*pm.AMU["F"]
            self.position = np.array([ [0,0,0],
                                      [0,0,1.554],
                                      [0,1.554,0],
                                      [1.554,0,0],
                                      [0,-1.554,0],
                                      [-1.554,0,0],
                                      [0,0,-1.554]])
            self.center = [0,0,0]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        elif self.name.upper() == "H2SO4":
            self.size = 7
            # Basic Parameters
            self.atom_names = ["S", "O", "O", "O", "O", "H", "H"]
            self.col_sphere = pm.Collision_Sphere["H2SO4"]
            self.mass = pm.AMU["S"] + 4*pm.AMU["O"] + 2*pm.AMU["H"]
            self.position = np.array([[0,0,0.1534],
                                      [0,1.2438, 0.8192	],
                                      [0,-1.2438, 0.8192],
                                      [ 1.2193	,0.0242,-0.8376],
                                      [-1.2193,-0.0242,-0.8376],
                                      [1.4709,-0.8641,-1.08],
                                      [-1.4709,0.8641, -1.0800]])
	
            self.center = [0,0,0]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        elif self.name.upper() == "CUBANE":
            self.size = 16
            self.atom_names = ["C"	,"C","C","C","C","C","C","C","H","H","H","H","H","H","H","H"]
            self.col_sphere = pm.Collision_Sphere["H2SO4"]
            self.position = np.array([[0.7854,0.7854,0.7854],
                                      [-0.7854,0.7854,0.7854],
                                      [0.7854,0.7854,-0.7854],
                                      [-0.7854,0.7854,-0.7854],
                                      [0.7854,-0.7854,0.7854],
                                      [-0.7854,-0.7854,0.7854],
                                      [0.7854,-0.7854,-0.7854],
                                      [-0.7854,-0.7854,-0.7854],
                                      [1.4188,1.4188,1.4188],
                                      [-1.4188,1.4188,1.4188],
                                      [1.4188,1.4188,-1.4188],
                                      [-1.4188,1.4188,-1.4188],
                                      [1.4188,-1.4188,1.4188],
                                      [-1.4188,-1.4188,1.4188],
                                      [1.4188,-1.4188,-1.4188],
                                      [-1.4188,-1.4188,-1.4188]])
            self.mass = 8*(pm.AMU["C"] + pm.AMU["H"])
            self.center = [0,0,0]
            self.com_frame = self.position.copy()
            self.previous_position = self.position.copy()
            self.type = "molecule"


        # Format of xyz file is assumed to be size, name, coordinates
        elif self.name.upper() == "CUSTOM":
            with open(custom_path, "r") as f:
                f1 = f.readlines()
                self.size = int(f1[0].split()[0])
                self.name = f1[1].split()[0]

                self.atom_names = []
                self.mass = 0

                for i in range(2, self.size+2):
                    Split = f1[i].split()
                    atom_name = Split[0]
                    self.atom_names.append(atom_name)
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
            self.type = "molecule"
            
        else:
            self.type = None

    def COM(self):
        R = self.position.copy()
        self.center = [0.0,0.0,0.0]
        A = 0
        for i in range(self.size):
            A += pm.AMU[self.atom_names[i]] * R[i] # m_i * r_i (Use AMU to get m_i and self.postion list comprehension to find r_i)
        self.center = A/self.mass


    def COM_FRAME(self):
        self.com_frame = 0
        R = deepcopy(self.position)
        R -= self.center * np.ones(R.shape)
        # for i in range(3):
        #     R[i] -= np.array(self.center)
            
        self.com_frame = R


    def translate(self,delta):
        delArr = np.array(delta)* np.ones(self.position.shape)
        self.position += delArr

    def xrotation(self, theta):

        theta = theta*math.pi/180 # converts theta to radians
        xvec = np.array([[1,0,0],
                        [0,np.cos(theta),np.sin(theta)],
                        [0,-np.sin(theta), np.cos(theta)]])

        R = np.dot(self.com_frame, xvec) # Find new matrix in COM Frame
        Delta = R - self.com_frame	   # Find changes in position from rotation


        self.position += Delta      			# Add difference to position matrix



    def yrotation(self, theta):
        theta = theta*np.pi/180.0 # converts theta to radians
        yvec = np.array([[np.cos(theta), 0 , np.sin(theta)],
                                         [0,1,0],
                                 [-np.sin(theta),0, np.cos(theta)]])

        R = np.dot(self.com_frame, yvec) # Find new matrix in COM Frame
        Delta = R - self.com_frame	   # Find changes in position from rotation

        self.position += Delta      			# Add difference to position matrix


    def zrotation(self, theta):
        theta = theta*math.pi/180 # converts theta Sto radians
        zvec = np.vstack([[np.cos(theta), -np.sin(theta), 0],
                                                 [np.sin(theta),np.cos(theta),0],
                                                 [0,0,1]])

        R = np.dot(self.com_frame, zvec) # Find new matrix in COM Frame
        Delta = R - self.com_frame	   # Find changes in position from rotation

        self.position += Delta      			# Add difference to position matrix



    def update(self):
        # self.previous_position = self.position.copy()
        self.COM()
        self.COM_FRAME()		
        
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

    # Copies previous position and then randomly moves the particle
    def random_move(self, dist=1, theta=90):
        if self.size == 1:
            self.translate([dist*(2*np.random.random()-1),  # Translate each randomly
                            dist*(2*np.random.random()-1),
                            dist*(2*np.random.random()-1)]) 		
            self.update()
        else:
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



class Create_System():
    
    def homogeneous(Npart, Mol_Class=Molecule, moleculeName=molName, create_xyz=True, save_path=filePath):
        M=[]
        for i in range(Npart):
            M.append(Mol_Class(moleculeName))  # Initialize the molecule	
            M[i].random_move(MAXT, MAXR) # Randomly move particle
    
            # This portion ensures molecules do not collide
            if i > 0: # Start checking after 1st particle
                j = 0 # Initialize j for while loop
                max_search = 100
                move_size = MAXT # Initial move_size is set by MAXT, can later be increased
                m = 0
                while j < i:
                    if calc.object_distance(M[i],M[j]) >= 2*COL_DIST: # Check if within collision distance
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
        if create_xyz:
            data.XYZ.create(M, path=save_path)
        return M    
    
    @classmethod 
    def extract_species(cls,species_list):
        name_list = []
        number_objects = []
        for i, name in enumerate(species_list):
            name_list.append(name)
            number_objects.append(species_list[name])
        return name_list, number_objects
    
    
    
    
    @classmethod    
    def heterogeneous(cls, species_list, Mol_Class=Molecule, Atom_Class= Atom, create_xyz=True, save_path=filePath):
        classifications = []
        names, number_obj = cls.extract_species(species_list)
        M = []
        collision_distance = []
        C = 0
        for k,name in enumerate(names):
            for i in range(number_obj[k]):
                if Mol_Class(name).type == "molecule":
                    classifications.append("molecule")
                    M.append(Molecule(name))
                elif Atom_Class(name).type == "atom":
                    classifications.append("atom")
                    M.append(Atom(name))
                else:
                    print("Error in Species List")
                    break 
        
                collision_distance.append(M[C].col_sphere)
                if M[C].size == 1:
                    M[C].random_move(MAXT)
                else:
                    M[C].random_move(MAXT, MAXR) # Randomly move particle
                
                if C > 0:
                    j = 0
                    max_search = MAX_SEARCH 
                    move_size = MAXT # Initial move_size is set by MAXT, can later be increased
                    m = 0
                    while j < C:
                        if calc.object_distance(M[C], M[j]) >= collision_distance[j] + collision_distance[C]: # Check if within collision distance
                            j += 1
                        else:
                            M[C].position = M[C].previous_position # Set equal to previous position
                            M[C].update() 
                            M[C].random_move(move_size, MAXR) # Randomly move again
                            j = 0
                            m += 1
                        if m > max_search:
                            move_size *= INCREASE_MOVE
                            m = 0                    
                C += 1
        if create_xyz:
            data.XYZ.create(M, path=save_path, homogeneous= False)
        
        return M    
    

