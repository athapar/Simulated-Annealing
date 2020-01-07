import pandas as pd
import numpy as np

# Can go to initialize
def Store_Data(M):
	X = pd.DataFrame(M[0].position[0])
	Y = pd.DataFrame(M[0].position[1])
	Z = pd.DataFrame(M[0].position[2])
	
	for i in range(1,len(M)):
		x = pd.DataFrame(M[i].position[0])
		y = pd.DataFrame(M[i].position[1])
		z = pd.DataFrame( M[i].position[2])
		
		X = pd.concat([X, x], axis = 0, ignore_index=True)
		Y = pd.concat([Y, y], axis = 0, ignore_index=True)
		Z = pd.concat([Z, z], axis = 0, ignore_index=True)
	return X, Y, Z

# Can go to initialize
def Append_Data(M, X_ar, Y_ar, Z_ar):
	Xi = pd.DataFrame(M[0].position[0])
	Yi = pd.DataFrame(M[0].position[1])
	Zi = pd.DataFrame(M[0].position[2])
	
	for i in range(1,len(M)):
		xi = pd.DataFrame(M[i].position[0])
		yi = pd.DataFrame(M[i].position[1])
		zi = pd.DataFrame(M[i].position[2])
		
		Xi = pd.concat([Xi, xi], axis = 0, ignore_index=True)
		Yi = pd.concat([Yi, yi], axis = 0, ignore_index=True)
		Zi = pd.concat([Zi, zi], axis = 0, ignore_index=True)
	
	X_ar = pd.concat([X_ar, Xi], axis = 1, ignore_index=True)
	Y_ar = pd.concat([Y_ar, Yi], axis = 1, ignore_index=True)
	Z_ar = pd.concat([Z_ar, Zi], axis = 1, ignore_index=True)
	
#	X_ar.index, Y_ar.index, Z_ar.index = range(1, len(X_ar.index)+1), range(1, len(Y_ar.index)+1), range(1, len(Z_ar.index)+1)
#	X_ar.columns, Y_ar.columns, Z_ar.columns = range(1, len(X_ar.columns)+1), range(1, len(Y_ar.columns)+1), range(1, len(Z_ar.columns)+1)
	return X_ar, Y_ar, Z_ar

def _to_array(mol_list):
	L = len(mol_list)
	for i in range(L):
		if i ==0:
			R = mol_list[i].position
		else:
			R = np.hstack([R, mol_list[i].position])
	
	return R

# Creates new xyz file with default name testprint			
def to_xyz(M, filename="testprint"):
	name = M[0].name
	f = open('Sim/%s.xyz' % filename,'w+')
	f.write("%i\n" % (M[0].size*len(M)))
	f.write("%s\n" % name)
	
	for i in range(len(M)):
		R = M[i].position
		X, Y, Z = R[0], R[1], R[2]
		for j in range(len(X)):
			f.write("%s \t %0.18g \t %0.18g \t %0.18g \n" % (M[i].mol_type[j],X[j], Y[j], Z[j]))
	f.close()
	return X, Y, Z

# Appends already existing xyz file. Default filename is testprint.
def append_xyz(M, filename="testprint"):
	name = M[0].name
	f = open("Sim/%s.xyz" % filename, "a")
	f.write("%i\n" % (3*len(M)))
	f.write("%s \n" % name)
	for i in range(len(M)):
		R = M[i].position
		X, Y, Z = R[0], R[1], R[2]
		for j in range(len(X)):
			f.write("%s \t %0.18g \t %0.18g \t %0.18g \n" % (M[i].mol_type[j],X[j], Y[j], Z[j]))
	f.close()
	return X, Y, Z
		
def create_df(M):
	R = []
	for i in range(len(M)):
		if i == 0:
			R = pd.DataFrame(M[i].position)
			T = pd.DataFrame(M[i].mol_type)
			R = pd.concat([T, R],axis=1,ignore_index=True)
		else:
			B = pd.DataFrame(M[i].position)
			T = pd.DataFrame(M[i].mol_type)
			B = pd.concat([T, B],axis=1,ignore_index=True)
			R = pd.concat([R,B], ignore_index = True)
	
	R.columns = ["Atom", "X", "Y", "Z"]
	R.index = range(1,len(R)+1)
	return R
	

