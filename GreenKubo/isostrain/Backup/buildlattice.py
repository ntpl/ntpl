import numpy as np
import itertools 

def build(n):
	lattice=np.array([[0,0,0]])
	basis=np.array([[0.0, 0.0, 0.0],
       [0.5, 0.5, 0.0],
       [0.5, 0.0, 0.5],
       [0.0, 0.5, 0.5]])
	for i in list(itertools.product(range(0,n),repeat=3)):
		x,y,z=i
		lattice=np.append(lattice,basis+[x,y,z],axis=0)
	lattice=np.delete(lattice,0,0)
	return lattice

	fcc = np.array([[0.0, 0.0, 0.0],
		[0.5, 0.5, 0.0],
		[0.5, 0.0, 0.5],
		[0.0, 0.5, 0.5]])

	


print build(2)

