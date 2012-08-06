import numpy as np
class Eigvec:

	def __init__(self,file_name):
		self._all=np.loadtxt(str(file_name),dtype=complex)

	def x(self,kindex=None,mode=None):
		if kindex == None and mode == None:
			return self._all[0::3,:]
		if mode != None:
			return self._all[0::3,mode]
		else:
			return self._all[len(self._all[0,:])*kindex:(len(self._all[0,:])*(kindex+1)-1):3,mode]

	def y(self,kindex=None,mode=None):
		if kindex == None and mode == None:
			return self._all[1::3,:]
		if mode != None:
			return self._all[1::3,mode]
		else:
			return self._all[len(self._all[0,:])*kindex+1:(len(self._all[0,:])*(kindex+1)):3,mode]

	def z(self,kindex=None,mode=None):
		if kindex == None and mode == None:
			return self._all[2::3,:]
		if mode != None:
			return self._all[2::3,mode]
		else:
			return self._all[len(self._all[0,:])*kindex+2:(len(self._all[0,:])*(kindex+1)+1):3,mode]

e=Eigvec('eigvec1.dat')
print e.y(mode=0)
