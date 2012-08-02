## ## ## strchg.py v.0.1
## ## ## This program simplifies and executes
## ## ## common text manipulation operations.
## ## ## Created: 07/29/2012 - KDP
## ## ## Last Edited: 07/31/2012 - KDP
import os

def sed(strings, orig, new):
	"""
	Simplifies and executes the sed terminal command

	ntplpy.strchg.sed(strings, orig, new)
	Parameters
	----------
		strings : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		orig : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		new : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
			An empty string activates the --in-place (-i)
			switch, which edits the original file without
			creating a new one.
	"""

	## Create concatenated string of command expressions
	commands = ''
	for key in strings:
		commands = commands +'-e \'s/' + str(key) + '/' + str(strings[key]) + '/g\' '
	
	## Execute the sed command for new file or --in-place
	if new == '':
		os.system('sed -i ' + commands + str(orig))
	else:
		os.system('sed ' + commands + str(orig) + ' > ' + str(new))
##### END SED


def grep(self):
	return None
