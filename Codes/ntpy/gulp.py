## ## ## gulp.py v.0.1
## ## ## This program executes simple gulp
## ## ## operations.
## ## ## Created: 10/31/2012 - KDP
## ## ## Last Edited: 10/31/2012 - KDP
import os
import numpy as np
import strchg as st
import param.lj as lj
import param.const as ct

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
	freq[:] = freq[:] * 2.0 * np.pi * ct.value('tocenti') * lj.value('tau') * ct.value('c') 
	# Remove the .freq file
	os.system('rm '+ tempName+ '.freq')

	return freq






