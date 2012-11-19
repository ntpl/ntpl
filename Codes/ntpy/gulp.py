## ## ## gulp.py v.0.1
## ## ## This program executes simple gulp
## ## ## operations.
## ## ## Created: 10/31/2012 - KDP
## ## ## Last Edited: 10/31/2012 - KDP
import os
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
	os.system('gulp <'+ tempName+ ' > output.gulp')
	# os.system('gulp '+ tempName+ ' output.gulp')
	# Extract the frequencies
	freq = np.zeros( (3 * numAtomsUC), dtype=float)
	freq = np.loadtxt(tempName+ '.freq', comments='--')
	freq[:] = freq[:] * 2.0 * np.pi * ct.value('c') * ct.value('tocenti') * lj.value('tau')
	# Remove the .freq file
	os.system('rm '+ tempName+ '.freq')

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
	os.system('gulp <'+ tempName+ ' > output.gulp')
	# Grep out eigenvectors
	st.grep(' 1 x', 3 * numAtomsUC, 'output.gulp', 'eigvec_grep.dat')
	xyzDict = dict({ 'x' : ''})
	st.sed(xyzDict, 'eigvec_grep.dat', 'eigvec2.dat')
	os.system('rm eigvec_grep.dat')
	xyzDict = dict({ 'y' : ''})
	st.sed(xyzDict, 'eigvec2.dat', 'eigvec3.dat')
	os.system('rm eigvec2.dat')
	xyzDict = dict({ 'z' : ''})
	st.sed(xyzDict, 'eigvec3.dat', 'eigvec4.dat')
	os.system('rm eigvec3.dat')
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
	
	os.system('rm eigvec4.dat')
	return eig

def vel(strchg, numAtomsUC, gulpName, tempName, kpt, deltaKpt=10e-5):
	"""
	Executes a gulp program with the string changes neccesary
	and returns the velocities.

	ntpy.gulp.vel(strchg, numAtomsUC, gulpName, tempName, deltaKpt=10e-5)
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
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName)
			# Change kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValMdk = freq(strchg, numAtomsUC, gulpName, tempName)
			vel[:, idim] = ((freqVal - freqValMdk) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		elif kpt[idim] == -0.5: # kpt at left boundary
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName)
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName)
			vel[:, idim] = ((freqValPdk - freqVal) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		elif kpt[idim] == 0.0: # kpt at gamma point
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName)
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName)
			vel[:, idim] = ((freqValPdk - freqVal) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)
		else: # kpt at gamma point
			freqVal = freq(strchg, numAtomsUC, gulpName, tempName)
			# Change kpt
			kpt[idim] = kpt[idim] + deltaKpt
			strchg = _kptstrchg(strchg, kpt)
			freqValPdk = freq(strchg, numAtomsUC, gulpName, tempName)
			kpt[idim] = kpt[idim] - (2.0 * deltaKpt)
			strchg = _kptstrchg(strchg, kpt)
			freqValMdk = freq(strchg, numAtomsUC, gulpName, tempName)
			vel[:, idim] = ((freqValPdk - freqValMdk) / deltaKpt) * 1e-10
			# Reset kpt
			kpt[idim] = kpt[idim] - deltaKpt
			strchg = _kptstrchg(strchg, kpt)

	return vel








