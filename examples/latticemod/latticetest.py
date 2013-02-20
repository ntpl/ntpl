## ## ## This script will provide examples of how to implement
## ## ## the ntpy lattice module

ntplPATH = '/home/kevin/projects/ntpl/Codes'

# First import the module
import sys
sys.path.append(ntplPATH) # Needed to link the ntplpy module
import ntpy.lattice as lt
# import numpy
import numpy as np


## ## Create single lammps file with single lattice.				## ##
## ## This is the most basic case possible.							## ##

#	First, create a lammps object that represents one file.
mylamp = lt.Lammps('test_1.txt') 
# mylamp = lt.lammps('filenamehere') 

# Next, create your lattice *NB: atom_mass, must be defined for lammps files.
mylat = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 4], atom_mass = [1])
# mylat = lt.block([lattice vector x, lat vec y, lat vec z], 'lattice type', [units in x, units in y, units in z], atom_mass = [atom mass here])

# Finally, build your lammps file using the block you just created.
mylamp.buildLammps([mylat])
## ## And you're done!!!

## ## We can also return the positions as a numpy array
npy_lattice = lt.buildNumpy([mylat])
print npy_lattice


## ## "blocks" are very modular, they can be called multiple times	## ##
## ## This creates the same lattice above but from two "blocks".	## ##

mylamp = lt.Lammps('test_2.txt')
# Notice this block is half the size in the z-direction
mylat = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 2], atom_mass = [1])
# Now pass the block twice to create two identical halves.
# The lattice will be concatenated in the z-direction automatically
# creating an identical lattice to the one above.
mylamp.buildLammps([mylat, mylat])
## ## End


## ## We can also build an xyz file for visualization				## ##
## ## by adding the atom_type parameter needed for xyz				## ##

myxyz = lt.Xyz('test_3.xyz')
# NB: atom_mass is not required for xyz files, therefore it is not included here
mylat = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 4], atom_type=['Ar'])
myxyz.buildXyz([mylat])
## ## End


## ## Another cool ability of the lattice module is					## ##
## ## multi-species lattices										## ##

mylamp = lt.Lammps('test_4.txt')
# Here we simply add the extra masses to out 'atom_mass' parameter to match the atoms
# in our basis vector. Therefore, our 'atom_mass' parameter must have the same length
# as the number of atoms in our basis vector
mylat = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 4], atom_mass = [1, 4, 1, 4])
mylamp.buildLammps([mylat])

# To see how the basis vector looks for a particular lat_type simply run:
lt.printBasis('fcc')
## ## End


## ## Let's do something more complicated.							## ##
## ## Here we create 4 files that alternate heave and light argon	## ##
## ## in increasing number. We also create xyz files to visualize	## ##
## ## our results from the same blocks.								## ##

hot_reservoir = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 2], atom_mass=[1], atom_type=['Ar'], bd_space=1.5632)
cold_reservoir = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 2], atom_mass=[1], atom_type=['Ar'])
# Note here we change the boundary spacing (bd_space), in this case we add (positive)
# one unit cell of spacing after the ORIGINAL location of the periodic boundary for
# the light and heavy argon blocks.
light_ar = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 4], atom_mass=[1], atom_type=['Ar'], bd_space=1.5632)
heavy_ar = lt.Block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 2], atom_mass=[4], atom_type=['Ar'], bd_space=1.5632)

for filenum in range(4):
	layers = [hot_reservoir]
	mylamp = lt.Lammps('alt_'+ str(filenum)+ '.lammps')
	myxyz = lt.Xyz('alt_'+ str(filenum)+ '.xyz')
	for lay in range(filenum + 1):
		layers.append(light_ar)
		layers.append(heavy_ar)
	layers.append(cold_reservoir) #end cap reservoir
	mylamp.buildLammps(layers)
	myxyz.buildXyz(layers)
## ## End

