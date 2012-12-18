## ## ## gulp.py v.0.1
## ## ## This program executes simple gulp
## ## ## operations.
## ## ## Created: 10/31/2012 - KDP
## ## ## Last Edited: 10/31/2012 - KDP
from os import system
import numpy as np
import strchg as st
import param.const as ct
import param.lj as lj

def freq(strchg, numAtomsUC, gulpName, tempName):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the frequencies.

	ntpy.gulp.freq(strchg, numAtomsUC, gulpName, tempName)
	Parameters
	----------
		strchg : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		numAtomsUC : int
			The number of atoms in the unit cell.
		gulpName : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		tempName : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
	"""
	# String change the template file
	st.sed(strchg, gulpName, tempName)
	# Execute the new gulp file
	system('gulp <'+ tempName+ ' > output.gulp')
	# system('gulp '+ tempName+ ' output.gulp')
	# Extract the frequencies
	freq = np.zeros( (3 * numAtomsUC), dtype=float)
	freq = np.loadtxt(tempName+ '.freq', comments='--')
	# Remove the .freq file
	system('rm '+ tempName+ '.freq')

	return freq

def eig(strchg, numAtomsUC, gulpName, tempName, kpt):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the eigenvectors.

	ntpy.gulp.eig(strchg, numAtomsUC, gulpName, tempName, kpt)
	Parameters
	----------
		strchg : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		numAtomsUC : int
			The number of atoms in the unit cell.
		gulpName : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		kpt : array of type float
			A numpy array of size 3 that has the three
			dimension k-point to be used.
		tempName : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
	"""
	# String change the template file
	st.sed(strchg, gulpName, tempName)
	# Execute the new gulp file
	system('gulp <'+ tempName+ ' > output.gulp')
	# Grep out eigenvectors
	st.grep(' 1 x', 3 * numAtomsUC, 'output.gulp', 'eigvec_grep.dat')
	xyzDict = dict({ 'x' : ''})
	st.sed(xyzDict, 'eigvec_grep.dat', 'eigvec2.dat')
	system('rm eigvec_grep.dat')
	xyzDict = dict({ 'y' : ''})
	st.sed(xyzDict, 'eigvec2.dat', 'eigvec3.dat')
	system('rm eigvec2.dat')
	xyzDict = dict({ 'z' : ''})
	st.sed(xyzDict, 'eigvec3.dat', 'eigvec4.dat')
	system('rm eigvec3.dat')
	# Read in eigvecs
	if kpt[0] == 0.0 and kpt[1] == 0.0 and kpt[2] == 0.0: 
		# If Gamma Point gulp only prints the real values
		dummy = np.loadtxt('eigvec4.dat', usecols=(1,2,3,4,5,6), comments='--')
		dummy = np.reshape(dummy, (3 * numAtomsUC, -1))
		eig = np.zeros( (3 * numAtomsUC, 3 * numAtomsUC), dtype=complex)
		for i in range(3 * numAtomsUC):
			eig[i,:] = dummy[i,:]
	else:
		eig = np.loadtxt('eigvec4.dat', usecols=(1,2,3,4,5,6), comments='--').view(complex)
		eig = np.reshape(eig, (3 * numAtomsUC, -1))
	
	system('rm eigvec4.dat')
	return eig

def vel(strchg, numAtomsUC, gulpName, tempName, kpt, freqConv=1.0, deltaKpt=10e-5):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the velocities.

	ntpy.gulp.vel(strchg, numAtomsUC, gulpName, tempName, freqConv=1.0, deltaKpt=10e-5)
	Parameters
	----------
		strchg : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		numAtomsUC : int
			The number of atoms in the unit cell.
		gulpName : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		tempName : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
		kpt : array of type float
			A numpy array of size 3 that has the three
			dimension k-point to be used.
		freqConv : float
			A float that can be used to convert the units of
			the freqency used. Defaults to 1.0 (no conversion).
		deltaKpt : float
			A float that determines what difference to use
			for the difference theorem. Defaults to 10e-5.
	"""
	def _kptstrchg(strchg, kpt):
		strchg['KPT'] = '{:10f} {:10f} {:10f}'.format(kpt[0], kpt[1], kpt[2])
		return strchg

	vel = np.zeros( (3 * numAtomsUC, 3), dtype=float)

	# For all three directions
	for idim in range(3):
		if kpt[idim] == 0.5: # kpt at right boundary
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			# Change kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValMdk = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			vel[:, idim] = ((freqVal - freqValMdk) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		elif kpt[idim] == -0.5: # kpt at left boundary
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			vel[:, idim] = ((freqValPdk - freqVal) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		elif kpt[idim] == 0.0: # kpt at gamma point
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			vel[:, idim] = ((freqValPdk - freqVal) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		else: # kpt at gamma point
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			kpt[idim] = kpt[idim] - (2.0 * deltaKpt)
			strchg = _kptstrchg(strchg, kpt)
			freqValMdk = freq(strchg, numAtomsUC, gulpName, tempName) * freqConv
			vel[:, idim] = ((freqValPdk - freqValMdk) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)

	return vel


def kpt(dim, conv=False):
	"""
	Creates an Lennard-Jones k-point list using either conventional or
	primative lattice vectors.

	ntpy.gulp.kpt(dim, conv = False)
	Parameters
	----------
		dim : array of type int
			An array that contains the number of unitcells in all three
			directions.
		conv : bool
			A boolean that describes whether to use the primative or
			conventional lattice vector. Defautls to conv = False, i.e. to use
			the primative lattice vector.
	"""

	# Total number of cells
	numdim = 1
	for i in range(len(dim)):
		numdim = numdim * dim[i]

	# Primative Lattice Vector
	a1 = np.array( [0.5, 0.5, 0], dtype=float )
	a2 = np.array( [0.5, 0, 0.5], dtype=float )
	a3 = np.array( [0, 0.5, 0.5], dtype=float )

	# Conventional Lattice Vector
	if conv is True: # Rewrite vectors
		a1 = np.array( [1.0, 0.0, 0.0], dtype=float )
		a2 = np.array( [0.0, 1.0, 0.0], dtype=float )
		a3 = np.array( [0.0, 0.0, 1.0], dtype=float )

	# Find reciprical lattice vector
	b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
	b2 = (np.cross(a1, a3)) / (np.dot(a2, np.cross(a1, a3)))
	b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))

	# The k-point matrix
	kpt = np.zeros( (numdim, 3) )

	## ## ## Create k-point list in integer format (i.e. missing 2pi/a)
	index = 0
	for ix in np.linspace(-(dim[0] / 2) + 1, dim[0] / 2, num = dim[0]):
		for iy in np.linspace(-(dim[1] / 2) + 1, dim[1] / 2, num = dim[1]):
			for iz in np.linspace(-(dim[2] / 2) + 1, dim[2] / 2, num = dim[2]):
				kpt[index, :] = (ix / dim[0]) * b1[:] + \
					(iy / dim[1]) * b2[:] + \
					(iz / dim[2]) * b3[:]
				index += 1	
	
	if conv is True: # No need to check against rec lat if conv
		return kpt

	# Build reciprical lattice points for comparison
	numdimrec = (dim[0]*2 + 1) * (dim[1]*2 + 1) * (dim[2]*2 + 1) 
	reclat = np.zeros( (numdimrec - 1, 3) )


	index = 0
	for ix in np.linspace(-dim[0], dim[0], num = dim[0]*2 + 1):
		for iy in np.linspace(-dim[1], dim[1], num = dim[1]*2 + 1):
			for iz in np.linspace(-dim[2], dim[2], num = dim[2]*2 + 1):
				reclat[index, :] = ix * b1[:] + iy * b2[:] + iz * b3[:]
				if all(reclat[index, :] != [0, 0, 0]): #make sure not first bz point
					index += 1

	# Now compare kpt to lattice points to map into first Brillouin Zone
	for ikpt in range(kpt[:, 0].size):
		while True:
			first = True # Bool to describe if kpt passed all reciprical latvecs

			for ireclat in range(reclat[:, 0].size):
				r_reclat2 = 0
				r_origin2 = 0

				for axis in range(3):
					r_origin = 0 - kpt[ikpt, axis]
					r_origin2 = r_origin2 + r_origin * r_origin
					r_reclat = reclat[ireclat, axis] - kpt[ikpt, axis]
					r_reclat2 = r_reclat2 + r_reclat * r_reclat

				if r_reclat2 < r_origin2:
					kpt[ikpt, :] = kpt[ikpt, :] - reclat[ireclat, :]
					first = False

			if first:
				break
	
	return kpt



