import transformations as tsf
import numpy as np
import math
f = open('Sim/testprint.xyz','w+')

# Need a way to transform molecules and extract the right amount of data from each array
# There is a confusion in the X,Y, and Z arrays showing all molecules and atoms
# The molecular array contains all the atoms in each molecule
# There is some parsing issue and this is bleeding into the rotation algorithm
# Need to check that the rotations rotate entire molecules correctly and not atoms or too few atoms


# Ways to go about rotations: If using X, Y, and Z then need way to know which atoms to extract and create a molecular matrix 

X = np.array([0.7493682, 0.0, -0.7493682])
Y = np.array([0, 0, 0])
Z = np.array([ 0.4424329, -0.1653507, .4424329])


def to_xyz(M, filename="testprint"):
	f = open('Sim/%s.xyz' % filename,'w+')
	f.write("%i\n" % (3*pm.N))
	f.write("Water\n")
	
	for i in range(len(M)):
		R = M[i].position
		X, Y, Z = R[0], R[1], R[2]
		for j in range(len(X)):
			f.write("%s \t %0.18g \t %0.18g \t %0.18g \n" % (M[i].mol_type[j],X[j], Y[j], Z[j]))
	f.close()
	return 0


X = vec[0]
Y = vec[1]
Z = vec[2]

f.write('%d\n' % n_part)     # Write number of atoms
f.write('water\n') 

f.write("H\t %.18g \t %.18g \t %.18g\n" % (X[0],Y[0],Z[0]))  
f.write("O\t %.18g \t %.18g \t %.18g\n" % (X[1],Y[1],Z[1]))  
f.write("H\t %.18g \t %.18g \t %.18g\n" % (X[2],Y[2],Z[2]))  



# print(tsf.xrotation(90, 0, X, Y, Z))