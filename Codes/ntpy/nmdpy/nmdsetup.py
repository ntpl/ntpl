#Author samuel.huberman@gmail.com

import numpy as np

class Pos:

	def __init__(self,file_name=None):
		self._file_name=file_name
		self._all=np.loadtxt(str(file_name),skiprows=17)

	def natoms(self):
		return len(self._all[:,0])

	def x(self):
		return self._all[:,2]

	def y(self):
		return self._all[:,3]

	def z(self):
		return self._all[:,4]
	def lx(self):
#		return np.maximum(self._all[:,2])
		return self._all[-1,2]
	def ly(self):
#               return np.maximum(self._all[:,3]) 
		return self._all[-1,3]
	def lz(self):
#               return np.maximum(self._all[:,4]) 
		return self._all[-1,4]
class Vel:

	def __init__(self,file_name=None):
		self._file_name=file_name
		self._all=np.loadtxt(str(file_name),comments='--')

	def vxts(self):
		return self._all[:,0]

	def vyts(self):
		return self._all[:,1]

	def vzts(self):
		return self._all[:,2]
	def nvel(self):
		return len(self._all[:,0])
class Eig:

	def __init__(self,file_name=None):
		self._all=np.loadtxt(str(file_name),dtype=complex)

	def x(self,kindex=None,mode=None):
		if kindex == None and mode == None:
			return self._all[0::3,:]
		if mode != None and kindex == None:
			return self._all[0::3,mode]
		else:
			return self._all[len(self._all[0,:])*kindex:(len(self._all[0,:])*(kindex+1)-1):3,mode]

	def y(self,kindex=None,mode=None):
		if kindex == None and mode == None:
			return self._all[1::3,:]
		if mode != None and kindex == None:
			return self._all[1::3,mode]
		else:
			return self._all[len(self._all[0,:])*kindex+1:(len(self._all[0,:])*(kindex+1)):3,mode]

	def z(self,kindex=None,mode=None):
		if kindex == None and mode == None:
			return self._all[2::3,:]
		if mode != None and kindex == None:
			return self._all[2::3,mode]
		else:
			return self._all[len(self._all[0,:])*kindex+2:(len(self._all[0,:])*(kindex+1)+1):3,mode]

	def natomucell(self):
		return len(self._all[0,:])/3
