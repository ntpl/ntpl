## ## ## create_crystal_file.py v.0.1					## ## ##
## ## ## This program creates a lattice data file ## ## ##
## ## ## WARNING: Unless you are adding extra			## ## ##
## ## ## functionality there is no need to edit		## ## ##
## ## ## this program!!!													## ## ##
## ## ## Created: 06/17/12 - KDP									## ## ##
## ## ## Last Edited: 06/18/12 - KDP							## ## ##

import sys

## ## ## Obtain user parameters										## ## ##
print 'Hello! This program creates a custom lattice structure'
print 'First we need to get some parameters from you'

lat_const = float(raw_input('Please enter your lattice constant: ').strip())

## ## ## Editors note: This would be a good place to insert a multi-atom function

lat_type = raw_input('Please enter your lattice type (e.g. fcc, bcc, diamond): ').strip().lower()
if (lat_type != 'fcc') & (lat_type != 'bcc') & (lat_type != 'diamond'):
	sys.exit('ERROR: Only fcc, bcc, and diamond functionality currently exits')

numx = int(raw_input('Please enter your x-dim unit cell length: ').strip())

numy = int(raw_input('Please enter your y-dim unit cell length: ').strip())

numz = int(raw_input('Please enter your z-dim unit cell length: ').strip())

file_type = raw_input('Please enter desired file type (lammps, xyz): ').strip().lower()
if (file_type != 'lammps') & (file_type != 'xyz'):
	sys.exit('ERROR: Must be either \'lammps\' or \'xyz\'')

if file_type == 'xyz':
	atom_type = raw_input('Please enter the atom type: ').strip()

if file_type == 'lammps':
	atom_mass = float(raw_input('Please enter the atom mass: ').strip())

file_name = raw_input('Please enter the output file name: ').strip()


## ## ## Begin print sequence												## ## ##

print_file = open(file_name, 'w')

atom_num = 1

dim = [numx, numy, numz]

total_atom_num = str(4 * dim[0] * dim[1] * dim[2])

## Entering LAMMPS print ##
if file_type == 'lammps':
	## Write LAMMPS format Header lines ##
	print_file.write('\n' + total_atom_num + ' atoms\n')
	print_file.write('1 atom types\n\n')

	if lat_type == 'fcc':
		fcc = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.0],
			[0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5] ])

		x_max = str(0.5 * lat_const + (0.5 + dim[0] - 1) * lat_const)
		y_max = str(0.5 * lat_const + (0.5 + dim[1] - 1) * lat_const)
		z_max = str(0.5 * lat_const + (0.5 + dim[2] - 1) * lat_const)

		print_file.write('0.0 ' + x_max + ' xlo xhi\n')
		print_file.write('0.0 ' + y_max + ' ylo yhi\n')
		print_file.write('0.0 ' + z_max + ' zlo zhi\n\n')

		print_file.write('Masses\n\n')
		print_file.write('1 ' + str(atom_mass) + '\n\n')
		print_file.write('Atoms\n\n')
		
		for x_count in range(dim[0]):
			for y_count in range(dim[1]):
				for z_count in range(dim[2]):
					for b in range(4):
						print_file.write(str(atom_num) + '	' + '1')
						print_file.write('	' + str(fcc[b][0] * lat_const + x_count * lat_const))
						print_file.write('	' + str(fcc[b][1] * lat_const + y_count * lat_const))
						print_file.write('	' + str(fcc[b][2] * lat_const + z_count * lat_const) + '\n')
						atom_num += 1

	if lat_type == 'bcc':
		bcc = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.5] ])

		x_max = str(0.5 * lat_const + (0.5 + dim[0]) * lat_const)
		y_max = str(0.5 * lat_const + (0.5 + dim[1]) * lat_const)
		z_max = str(0.5 * lat_const + (0.5 + dim[2]) * lat_const)

		print_file.write('0.0 ' + x_max + ' xlo xhi\n')
		print_file.write('0.0 ' + y_max + ' ylo yhi\n')
		print_file.write('0.0 ' + z_max + ' zlo zhi\n\n')
		
		print_file.write('Masses\n\n')
		print_file.write('1 ' + str(atom_mass) + '\n\n')
		print_file.write('Atoms\n\n')
		
		for x_count in range(dim[0]):
			for y_count in range(dim[1]):
				for z_count in range(dim[2]):
					for b in range(2):
						print_file.write(str(atom_num) + '	' + '1')
						print_file.write('	' + str(bcc[b][0] * lat_const + x_count * lat_const))
						print_file.write('	' + str(bcc[b][1] * lat_const + y_count * lat_const))
						print_file.write('	' + str(bcc[b][2] * lat_const + z_count * lat_const) + '\n')
						atom_num += 1

	if lat_type == 'diamond':
		diamond = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.0],
			[0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5],
			[0.25, 0.25, 0.25],
			[0.75, 0.75, 0.25],
			[0.75, 0.25, 0.75],
			[0.25, 0.75, 0.75] ])

		x_max = str(0.75 * lat_const + (0.25 + dim[0]) * lat_const)
		y_max = str(0.75 * lat_const + (0.25 + dim[1]) * lat_const)
		z_max = str(0.75 * lat_const + (0.25 + dim[2]) * lat_const)

		print_file.write('0.0 ' + x_max + ' xlo xhi\n')
		print_file.write('0.0 ' + y_max + ' ylo yhi\n')
		print_file.write('0.0 ' + z_max + ' zlo zhi\n\n')
		
		print_file.write('Masses\n\n')
		print_file.write('1 ' + str(atom_mass) + '\n\n')
		print_file.write('Atoms\n\n')
		
		for x_count in range(dim[0]):
			for y_count in range(dim[1]):
				for z_count in range(dim[2]):
					for b in range(8):
						print_file.write(str(atom_num) + '	' + '1')
						print_file.write('	' + str(diamond[b][0] * lat_const + x_count * lat_const))
						print_file.write('	' + str(diamond[b][1] * lat_const + y_count * lat_const))
						print_file.write('	' + str(diamond[b][2] * lat_const + z_count * lat_const) + '\n')

## Entering xyz print ##
if file_type == 'xyz':
	## Write xyz format Header lines ##
	print_file.write(total_atom_num + '\n')
	print_file.write(file_name + '\n')

	if lat_type == 'fcc':
		fcc = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.0],
			[0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5] ])
		
		x_max = str(0.5 * lat_const + (0.5 + dim[0]) * lat_const)
		y_max = str(0.5 * lat_const + (0.5 + dim[1]) * lat_const)
		z_max = str(0.5 * lat_const + (0.5 + dim[2]) * lat_const)

		for x_count in range(dim[0]):
			for y_count in range(dim[1]):
				for z_count in range(dim[2]):
					for b in range(4):
						print_file.write(atom_type + '	' + str(atom_num) + '	' + str(atom_mass))
						print_file.write('	' + str(fcc[b][0] * lat_const + x_count * lat_const))
						print_file.write('	' + str(fcc[b][1] * lat_const + y_count * lat_const))
						print_file.write('	' + str(fcc[b][2] * lat_const + z_count * lat_const) + '\n')
						atom_num += 1

	if lat_type == 'bcc':
		bcc = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.5] ])

		x_max = str(0.5 * lat_const + (0.5 + dim[0]) * lat_const)
		y_max = str(0.5 * lat_const + (0.5 + dim[1]) * lat_const)
		z_max = str(0.5 * lat_const + (0.5 + dim[2]) * lat_const)

		for x_count in range(dim[0]):
			for y_count in range(dim[1]):
				for z_count in range(dim[2]):
					for b in range(2):
						print_file.write(str(atom_num) + '	' + str(atom_mass))
						print_file.write('	' + str(bcc[b][0] * lat_const + x_count * lat_const))
						print_file.write('	' + str(bcc[b][1] * lat_const + y_count * lat_const))
						print_file.write('	' + str(bcc[b][2] * lat_const + z_count * lat_const) + '\n')
						atom_num += 1

	if lat_type == 'diamond':
		diamond = ([ [0.0, 0.0, 0.0],
			[0.5, 0.5, 0.0],
			[0.5, 0.0, 0.5],
			[0.0, 0.5, 0.5],
			[0.25, 0.25, 0.25],
			[0.75, 0.75, 0.25],
			[0.75, 0.25, 0.75],
			[0.25, 0.75, 0.75] ])

		x_max = str(0.75 * lat_const + (0.25 + dim[0]) * lat_const)
		y_max = str(0.75 * lat_const + (0.25 + dim[1]) * lat_const)
		z_max = str(0.75 * lat_const + (0.25 + dim[2]) * lat_const)

		for x_count in range(dim[0]):
			for y_count in range(dim[1]):
				for z_count in range(dim[2]):
					for b in range(8):
						print_file.write(str(atom_num) + '	' + str(atom_mass))
						print_file.write('	' + str(diamond[b][0] * lat_const + x_count * lat_const))
						print_file.write('	' + str(diamond[b][1] * lat_const + y_count * lat_const))
						print_file.write('	' + str(diamond[b][2] * lat_const + z_count * lat_const) + '\n')
##	##	##

print('Total number of atoms ' + str(total_atom_num))
print('Box size: [ ' + x_max + ', ' + y_max + ', ' + z_max + ' ]')
print('Program complete: Thank you for shopping at Wal-Mart')
