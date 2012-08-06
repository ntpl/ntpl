## ## ## param.py v.0.1
## ## ## This program returns common parameters
## ## ## Created: 08/05/2012 - KDP
## ## ## Last Edited: 08/05/2012 - KDP

def lj(string):
	"""
	Returns common parameters for lj argon.

	ntpy.param.lj(string)
	Parameters
	----------
		string : str
			A string that corresponds to a lennard jones argon parameter.
	"""

	ljparams = dict({
					'lat0' : 5.269,		# Lattice constant for argon at 0 K in angstroms
					'lat20' : 5.315,	# Lattice constant for argon at 20 K in angstroms
					'lat35' : 5.355,	# Lattice constant for argon at 35 K in angstroms
					'lat50' : 5.401,	# Lattice constant for argon at 50 K in angstroms
					'lat65' : 5.455,	# Lattice constant for argon at 65 K in angstroms
					'lat80' : 5.527,	# Lattice constant for argon at 80 K in angstroms
					'epsilon' : 1.67e-21,	# Epsilon constant for argon in joules
					'sigma' : 3.40e-10,	# Sigma constant for argon in meters
					'mass' : 6.63e-26,	# Mass constant for argon in kilograms
					'tau' : 2.14e-12	# Tau constant for argon in seconds
					})	

	try:
		return ljparams[string]
	except KeyError, e:
		print "KeyError: %s is not a valid key for ntpy.param.lj()." % e
		raise	
##### END LJ

def const(string):
	"""
	Returns common constant parameters.

	ntpy.param.const(string)
	Parameters
	----------
		string : str
			A string that corresponds to a constant parameter.
	"""

	constparams = dict({
					'kb' : 1.3806e-23,	# Boltzmann's constant
					'hbar' : 1.054e-34,	# Planck's constant
					'topeta' : 1e15,	# To peta-
					'totera' : 1e12,	# To tera-
					'togiga' : 1e9,	# To giga-
					'tomega' : 1e6,	# To mega-
					'tokilo' : 1e3,	# To kilo-
					'tocenti' : 1e-2,	# To centi-
					'tomilli' : 1e-3,	# To milli-
					'tomicro' : 1e-6,	# To micro-
					'tonano' : 1e-9,	# To nano-
					'topico' : 1e-12,	# To pico-
					'tofemto' : 1e-15,	# To femto-
					})

	try:
		return constparams[string]
	except KeyError, e:
		print "KeyError: %s is not a valid key for ntpy.param.const()." % e
		raise
	
##### END LJ

