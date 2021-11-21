""" This file collects data from molecular arrays and stores them into either xyz files to dataframes
    This file also contains function used in displaying potential energies and other data formatting needs"""

import pandas as pd
import numpy as np
import parameters as pm
import os
import matplotlib.pyplot as plt

defaultSavePath = os.path.join(os.getcwd(),"Sim", "testPrint")


class XYZ():
    """ Used to create and append xyz files """	
    def create(M, path=defaultSavePath, homogeneous = True):
        if homogeneous: 
            name = M[0].name
            size = sum([mol.size for mol in M])
            f = open(path+".xyz",'w+')
            f.write(str(int(size))+"\n")
            f.write(str(name)+"\n")
        else:
            size = sum([mol.size for mol in M]) # sum up the sizes to get total # atoms
            comment = "Mixed Species Sim"
            f = open(path+".xyz",'w+')
            f.write(str(int(size))+"\n")
            f.write(comment+"\n")
            
        for obj in M:
            if obj.size == 1:
                f.write("{} \t {:.6f} \t {:.6f} \t {:.6f} \n".format(obj.atom_names[0], obj.position[0], obj.position[1], obj.position[2]))
            else:
                for j,coords in enumerate(obj.position):
                    f.write("{} \t {:.6f} \t {:.6f} \t {:.6f} \n".format(obj.atom_names[j], coords[0], coords[1], coords[2]))
        f.close()
        

    # Appends already existing xyz file. Default filename is testprint.
    def append(M, path=defaultSavePath, homogeneous = True):
        if homogeneous:  
            name = M[0].name
            size = sum([mol.size for mol in M])
            f = open(path+".xyz", "a")
            f.write(str(int(size))+"\n")
            f.write(str(name)+"\n")
        else:
            size = sum([mol.size for mol in M]) # sum up the sizes to get total # atoms
            comment = "Mixed Species System"
            f = open(path+".xyz", "a")
            f.write(str(int(size))+"\n")
            f.write(comment+"\n")
            
        for obj in M:
            if obj.size == 1:
                f.write("{} \t {:.6f} \t {:.6f} \t {:.6f} \n".format(obj.atom_names[0], obj.position[0], obj.position[1], obj.position[2]))
            else:
                for iter1,coords in enumerate(obj.position):
                    f.write("{} \t {:.6f} \t {:.6f} \t {:.6f} \n".format(obj.atom_names[iter1], coords[0], coords[1], coords[2]))
        
        f.close()
    	

    
def create_df(M):
    """ Creates a dataframe for a given system """	
    R = []
    for i,obj in enumerate(M):
        if i == 0:
            if obj.type == "molecule":
                R = pd.DataFrame(obj.position)
                T = pd.DataFrame(obj.atom_names)
                R = pd.concat([T, R],axis=1,ignore_index=True)
            else:
                B = list(obj.position)
                B.insert(0,obj.name)
                R = pd.DataFrame(B).T
                
        else:
            if obj.type == "molecule":
                B = pd.DataFrame(M[i].position)
                T = pd.DataFrame(M[i].atom_names)
                B = pd.concat([T, B],axis=1,ignore_index=True)
                R = pd.concat([R,B], ignore_index = True)
            else:
                
                B = list(obj.position)
                B.insert(0,obj.name)
                B = pd.DataFrame(B).T
                R = pd.concat([R, B], ignore_index = True)
                
	
    R.columns = ["Atom", "X", "Y", "Z"]
    R.index = range(1,len(R)+1)
    return R
	
def create_dictionary(loc ="Spec_List.txt"):
    """ Creates dictionary from a text file: notation and example given in  Spec_List.txt"""	
    with open(loc,mode='r') as f:
        f1 = f.readlines()[1:] # 1st line (comment) omitted
        x = [line.replace(",","").split() for line in f1] # List comprehension to store each line within list
        x = list(filter(None,x)) # Removes blank lines
        dictionary = {element[0]: int(element[1]) for element in x} # Dictionary comprehension converts list to dictionary
    return dictionary


def createPaths(PATHS):
    """ Checks if paths exist or creates them if needed for saving files"""
    for path in PATHS:
        try:
            os.mkdir(path)
        except:
            pass
    return None


def create_storage_dict(U_list, names_list):
    try:
        len(U_list)
        return {name: u for (name, u) in zip(names_list, U_list)}
    except:
        return {names_list: U_list}
        
    
def names(M):
    names_list = []
    for obj in M:
        names_list.append(obj.name)
    return names_list 


def functionNames(funcs):
    names = []
    for func in funcs:
        names.append(func.__name__)
    return names


