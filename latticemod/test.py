## ## ## This script will provide examples of how to implement
## ## ## the lattice module

# First import the module
import lattice as lt

## ## Create single lammps file with single lattice.
## ## This is the most basic case possible.

#	First, create a lammps object that represents one file.
mylamp = lt.lammps('test_1.txt') 
# mylamp = lt.lammps('filenamehere') 

# Next, create your lattice *NB: atom_mass, must be defined for lammps files.
mylat = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 4], atom_mass = [1])
# mylat = lt.block([lattice vector x, lat vec y, lat vec z], 'lattice type', [units in x, units in y, units in z], atom_mass = [atom mass here])

# Finally, build your lammps file using the block you just created.
mylamp.buildlammps([mylat])
## ## And you're done!!!


## ## "blocks" are very modular, they can be called multiple times
## ## This creates the same lattice above but from two "blocks".

mylamp = lt.lammps('test_2.txt')
# Notice this block is half the size in the z-direction
mylat = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 2], atom_mass = [1])
# Now pass the block twice to create two identical halves.
# The lattice will be concatenated in the z-direction automatically
# creating an identical lattice to the one above.
mylamp.buildlammps([mylat, mylat])
## ## End


## ## We can also build an xyz file for visualization
## ## by adding the atom_type parameter needed for xyz

myxyz = lt.xyz('test_3.xyz')
# NB: atom_mass is not required for xyz files, therefore it is not included here
mylat = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 4], atom_type=['Ar'])
myxyz.buildxyz([mylat])
## ## End


## ## Another cool ability of the lattice module is multi-species lattices

mylamp = lt.xyz('test_4.xyz')
# Here we simply add the extra masses to out 'atom_mass' parameter to match the atoms
# in our basis vector. Therefore, our 'atom_mass' parameter must have the same length
# as the number of atoms in our basis vector
mylat = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [4, 4, 4], atom_mass = [1, 4, 1, 4])
mylamp.buildlamp([mylat])

# To see how the basis vector looks for a particular lat_type simply run:
lt.print_basis('fcc')
## ## End

## ## Let's do something more complicated
## ## Here we create 4 files that alternate heave and light argon in increasing
## ## number. We also create xyz files to visualize our results from the same
## ## blocks.

hot_reservoir = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 2], atom_mass=[1], atom_type=['Ar'], bd_space=1.5632)
cold_reservoir = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 2], atom_mass=[1], atom_type=['Ar'])
# Note here we change the boundary spacing (bd_space), in this case we add (positive)
# one unit cell of spacing after the ORIGINAL location of the periodic boundary for
# the light and heavy argon blocks.
light_ar = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 4], atom_mass=[1], atom_type=['Ar'], bd_space=1.5632)
heavy_ar = lt.block([1.5632, 1.5632, 1.5632], 'fcc', [8, 8, 2], atom_mass=[4], atom_type=['Ar'], bd_space=1.5632)

for filenum in range(4):
	layers = [hot_reservoir]
	mylamp = lt.lammps('alt_'+ str(filenum)+ '.lammps')
	myxyz = lt.xyz('alt_'+ str(filenum)+ '.xyz')
	for lay in range(filenum + 1):
		layers.append(light_ar)
		layers.append(heavy_ar)
	layers.append(cold_reservoir) #end cap reservoir
	mylamp.buildlammps(layers)
	myxyz.buildxyz(layers)
## ## End

