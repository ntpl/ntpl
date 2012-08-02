## ## ## lattice.py v.0.1
## ## ## This program creates lattice data files
## ## ## Created: 06/21/2012 - KDP
## ## ## Last Edited: 06/25/2012 - KDP

from sys import exit

## ## ## lammps class
class Lammps:
	"""
	Creates an object that represents a single lammps file.

	lattice.Lammps(file_name)
	Parameters
	----------
		file_name : str
			The name of the lammps file
	"""
	def __init__(self, file_name):
		self.file_name = file_name
	#-- END __init__ --#

	def buildLammps(self, blocks):
		"""
		This creates and writes the lammps file based on the blocks provided

		lattice.Lammps.buildLammps(blocks)
		Parameters
		----------
			blocks : list of type Block
				The blocks used to build the lammps file. Blocks
				should be put in order of origin to z-direction
				max.
		"""
		print_file = open(self.file_name, 'w')

		def _outPos(lat_vector, basis, atom_mass, dim, bd_space, origin, masses, atom_num):
			for x_count in range(dim[0]):
				for y_count in range(dim[1]):
					for z_count in range(dim[2]):
						for b in range(len(basis)):
							# Print atom number
							print_file.write(str(atom_num) + '\t')
							# Print atom type
							if len(atom_mass) == 1: #singular atom type
								print_file.write(str(masses[atom_mass[0]]) + '\t')
							else: 
								print_file.write(str(masses[atom_mass[b]]) + '\t')
							# Print positions	
							print_file.write(str(origin[0] + (basis[b][0] * lat_vector[0] + x_count * lat_vector[0])))
							print_file.write('\t')
							print_file.write(str(origin[1] + (basis[b][1] * lat_vector[1] + y_count * lat_vector[1])))
							print_file.write('\t')
							print_file.write(str(origin[2] + (basis[b][2] * lat_vector[2] + z_count * lat_vector[2])))
							print_file.write('\n')
							atom_num += 1

			origin[2] = (lat_vector[2] * dim[2]) + (origin[2] + bd_space) #take max length and add bdspace
			
			return origin, atom_num
		#-- END outPos --#


		## ## Preperation sequence
		origin = [0, 0, 0] #inital origin
		masses = dict() #list to hold the index of masses for all blocks
		mass_index = 1 #counter to create new atom/mass types
		atom_num = 1 #atom num/index
		total_atom_num = 0 #initalize total atom num
		x_max = 0 #initalize x_max
		y_max = 0 #initalize y_max
		z_max = 0 #initalize z_max
		
		## Parameter check
		for i in range(len(blocks)):
			if blocks[i].atom_mass == None:
				exit("ERROR: atom_mass parameter for "+ str(blocks[i])+ " block does not exist.")
			if len(blocks[i].atom_mass) > 1:
				if len(blocks[i].atom_mass) != len(blocks[i].basis):
					exit("ERROR: When declaring multiple species, the \"atom_mass\" parameter must be\n"+
							"the same length as the number of atoms in the basis vector. For more info\n"+
							"check the documentation (also lammps.print_basis).")

		## Loop over all Blocks
		for i in range(len(blocks)):
			basis_atom_num = len(blocks[i].basis) #how many atoms are in the basis
			total_atom_num = total_atom_num + (basis_atom_num * blocks[i].dim[0] * blocks[i].dim[1] * blocks[i].dim[2]) #total num of atoms
			
			# Find box boundaries
			x_temp = (blocks[i].lat_vector[0] * blocks[i].dim[0]) #finds largest x coordinate
			if x_temp > x_max:
				x_max = x_temp
			y_temp = (blocks[i].lat_vector[1] * blocks[i].dim[1]) #finds largest y coordinate
			if y_temp > y_max:
				y_max = y_temp
			z_max = z_max + (blocks[i].lat_vector[2] * blocks[i].dim[2]) + blocks[i].bd_space #lattices are stacked in zed

			# Create dict of mass:index values
			for mass in range(len(blocks[i].atom_mass)):
				if blocks[i].atom_mass[mass] not in masses:
					masses[blocks[i].atom_mass[mass]] = mass_index
					mass_index += 1

		## ## Print sequence	
		print_file.write('\n' + str(total_atom_num) + ' atoms\n') #print total num of atoms
		print_file.write(str(len(masses))+ ' atom types\n\n') #print num of atom types

		print_file.write('0.0 ' + str(x_max) + ' xlo xhi\n') #x-box size
		print_file.write('0.0 ' + str(y_max) + ' ylo yhi\n') #y-box size
		print_file.write('0.0 ' + str(z_max) + ' zlo zhi\n\n') #z-box size

		# Print masses
		print_file.write('Masses\n\n')
		for i, m in zip(range(len(masses)), iter(masses)):
			print_file.write(str(i + 1)+ '\t'+ str(m)+ '\n')
		print_file.write('\nAtoms\n\n')
	
		for i in range(len(blocks)):
			origin, atom_num = _outPos(blocks[i].lat_vector, blocks[i].basis, blocks[i].atom_mass, blocks[i].dim, blocks[i].bd_space, origin, masses, atom_num) #returns new origin

		print('Lammps file complete: Atoms '+ str(total_atom_num)+ '; Box size ['+ str(x_max)+ ', '+ str(y_max)+ ', '+ str(z_max)+ ']\n')
	#-- END buildLammps --#
#-- END Lammps --#

## ## ## xyz class
class Xyz:
	"""
	Creates an object that represents a single xyz file.

	lattice.Xyz(file_name)
	Parameters
	----------
		file_name : str
			The name of the xyz file
	"""
	def __init__(self, file_name):
		self.file_name = file_name
	#-- END __init__ --#

	def buildXyz(self, blocks):
		"""
		This creates and writes the xyz file based on the blocks provided

		lattice.Xyz.buildXyz(blocks)
		Parameters
		----------
			blocks : list of type block
				The blocks used to build the lammps file. Blocks
				should be put in order of origin to z-direction
				max.
		"""
		print_file = open(self.file_name, 'w')

		def _outPos(lat_vector, basis, dim, atom_type, bd_space, origin):
			for x_count in range(dim[0]):
				for y_count in range(dim[1]):
					for z_count in range(dim[2]):
						for b in range(len(basis)):
							# Print atom type
							if len(atom_type) == 1: #singular atom type
								print_file.write(str(atom_type[0]) + '\t')
							else:
								print_file.write(str(atom_type[b]) + '\t')
							# Print positions	
							print_file.write(str(origin[0] + (basis[b][0] * lat_vector[0] + x_count * lat_vector[0])))
							print_file.write('\t')
							print_file.write(str(origin[1] + (basis[b][1] * lat_vector[1] + y_count * lat_vector[1])))
							print_file.write('\t')
							print_file.write(str(origin[2] + (basis[b][2] * lat_vector[2] + z_count * lat_vector[2])))
							print_file.write('\n')

			origin[2] = (lat_vector[2] * dim[2]) + (origin[2] + bd_space) #take max length and add bdspace
			
			return origin
		#-- END outPos --#

		## ## Preperation sequence
		origin = [0, 0, 0] #inital origin
		total_atom_num = 0 #initalize total atom num
		x_max = 0 #initalize x_max
		y_max = 0 #initalize y_max
		z_max = 0 #initalize z_max

		## Parameter check
		for i in range(len(blocks)):
			if blocks[i].atom_type == None:
				exit("ERROR: atom_type parameter for "+ str(blocks[i])+ " block does not exist.")
			if len(blocks[i].atom_type) > 1:
				if len(blocks[i].atom_type) != len(blocks[i].basis):
					exit("ERROR: When declaring multiple species, the \"atom_type\" parameter must be\n"+
							"the same length as the number of atoms in the basis vector. For more info\n"+
							"check the documentation (also lattice.print_basis).")

		## Loop over all Blocks
		for i in range(len(blocks)):
			basis_atom_num = len(blocks[i].basis) #how many atoms are in the basis
			total_atom_num = total_atom_num + (basis_atom_num * blocks[i].dim[0] * blocks[i].dim[1] * blocks[i].dim[2]) #total num of atoms
			
		## ## Print sequence	
		print_file.write(str(total_atom_num) + ' atoms\n') #print total num of atoms
		print_file.write(self.file_name + '\n') #print file_name in comment line

		for i in range(len(blocks)):
			origin = _outPos(blocks[i].lat_vector, blocks[i].basis, blocks[i].dim, blocks[i].atom_type, blocks[i].bd_space, origin) #returns new origin

		print('XYZ file complete: Atoms '+ str(total_atom_num)+ '\n')
	#-- END buildXyz --#
#-- END Xyz --#


## ## ## Block class 
class Block:
	"""
	Defines a singular lattice that can linearly interfaced with other lattices

	lattice.Block(lat_vector, lat_type, dim, bd_space=0, atom_mass=None, atom_type=None)
	Parameters
	----------
		lat_vector : list of type float
			The lattice constants for the x, y, and z directions
		lat_type : str
			The type of lattice to be constructed. Currently only
			'sc', 'fcc', 'bcc', and 'diamond' are supported
		dim : list of type int
			The dimensions of the lattice in the x, y, and z directions.
			Dimensions are in units of lattice units.
		bd_space : float, optional
			The amount of space added at the terminating z boundary of
			the lattice. The space added is relative to the standard
			periodic foundary of the lattice. A negative value reduces
			the amount of space at the boundary.
		atom_mass : list of type float, optional
			The atom masses corresponding to the atoms in the basis vector.
			If only one mass is used then only one mass need be supplied.
			For a listing of the basis vector for each lattice type used
			by this module, use the lattice function "basis_print".
			Required for lammps files.
		atom_type : list of type str, optional
			The atom types of the lattice. Required for xyz files.
	"""
	def __init__(self, lat_vector, lat_type, dim, bd_space=0, atom_mass=None, atom_type=None):
		self.lat_vector = lat_vector
		self.lat_type = lat_type
		self.dim = dim
		self.bd_space = bd_space
		self.atom_mass = atom_mass
		self.atom_type = atom_type
		self.basis = []
		if lat_type=='sc' or lat_type=='cP' or lat_type=='cP1' or lat_type=='simple cubic' or lat_type=='primitive cubic':
			self.basis = ([ [0.0, 0.0, 0.0] ])
		elif lat_type=='bcc' or lat_type=='cI' or lat_type=='cI2' or lat_type=='body centered cubic':
			self.basis = ([ [0.0, 0.0, 0.0],
				[0.5, 0.5, 0.5] ])
		elif lat_type=='fcc' or lat_type=='cF' or lat_type=='cF4' or lat_type=='face centered cubic':
			self.basis = ([ [0.0, 0.0, 0.0],
				[0.5, 0.5, 0.0],
				[0.5, 0.0, 0.5],
				[0.0, 0.5, 0.5] ])
		elif lat_type=='diamond' or lat_type=='cF8':
			self.basis = ([ [0.0, 0.0, 0.0],
				[0.5, 0.5, 0.0],
				[0.5, 0.0, 0.5],
				[0.0, 0.5, 0.5],
				[0.25, 0.25, 0.25],
				[0.75, 0.75, 0.25],
				[0.75, 0.25, 0.75],
				[0.25, 0.75, 0.75] ])
		else:
			exit("ERROR: Only types \'sc\', \'fcc\', \'bcc\', \'diamond\' are supported")

	#-- END __init__ --#
#-- END Block --#

def printBasis(lat_type):
	"""
	Prints the basis vectors defined in this module for reference

	lattice.printBasis(lat_type)
	Parameters
	----------
		lat_type : str
			The type of lattice to be constructed. Currently only
			'sc', 'fcc', 'bcc', and 'diamond' are supported
	"""
	basis = []
	if lat_type=='sc' or lat_type=='cP' or lat_type=='cP1' or lat_type=='simple cubic' or lat_type=='primitive cubic':
		basis = ([ [0.0, 0.0, 0.0] ])
	elif lat_type=='bcc' or lat_type=='cI' or lat_type=='cI2' or lat_type=='body centered cubic':
		basis = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.5] ])
	elif lat_type=='fcc' or lat_type=='cF' or lat_type=='cF4' or lat_type=='face centered cubic':
		basis = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.0],
			[0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5] ])
	elif lat_type=='diamond' or lat_type=='cF8':
		basis = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.0],
			[0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5],
			[0.25, 0.25, 0.25],
			[0.75, 0.75, 0.25],
			[0.75, 0.25, 0.75],
			[0.25, 0.75, 0.75] ])
	else:
		exit("ERROR: Only types \'sc\', \'fcc\', \'bcc\', \'diamond\' are supported")

	for b in range(len(basis)):
		print 'Atom '+ str(b + 1)+ '\t'+ str(basis[b])
#-- END print_basis --#
	
	
