## ## ## Kevin Parrish 5/16/12 LAMMPS output editor ## ## ##

for a in range(5):
	orig = open('/home/kevin/jobs/GreenKubo/ljcrystal/relaxed_15K_' + str(a + 1) + '.pos', 'r')
	new = open('/home/kevin/jobs/GreenKubo/ljcrystal/in_structure_' + str(a + 1) + '.pos', 'w')

	## Read in Header lines ##
	junk = orig.readline()
	junk = orig.readline()
	junk = orig.readline()
	NumAtomsStr = orig.readline()
	junk = orig.readline()
	BoxBoundX = orig.readline()
	BoxBoundY = orig.readline()
	BoxBoundZ = orig.readline()
	junk = orig.readline()
	## ## ## ## ## ## ## ## ## ##

	## Write LAMMPS format Header lines ##
	new.write('\n' + NumAtomsStr[:-1] + ' atoms\n')
	new.write('1 atom types\n\n')
	
	new.write(BoxBoundX[:-1] + ' xlo xhi\n')
	new.write(BoxBoundY[:-1] + ' ylo yhi\n')
	new.write(BoxBoundZ[:-1] + ' zlo zhi\n\n')
	
	new.write('Masses\n\n1 1.0\n\nAtoms\n\n')
	## ## ## ## ## ## ## ## ## ##
	
	NumAtomsInt = int(float(NumAtomsStr)) #convert str to int
	
	## Write LAMMPS format Position lines ##
	for b in range(NumAtomsInt):
		curr = orig.readline()
		new.write(str(b + 1) + ' 1 ' + curr)
	## ## ## ## ## ## ## ## ## ##
	
	orig.close
	new.close
	
