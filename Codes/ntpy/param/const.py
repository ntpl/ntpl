## ## ## const.py v.0.1
## ## ## This program returns common parameters
## ## ## Created: 08/07/2012 - KDP

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

