#Author samuel.huberman@gmail.com

import numpy as np

class Nmd:

	def __init__(self,vel,pos,eig,kpt):
			self._vel=vel
			self._pos=pos
			self._eig=eig
			self._kpt=kpt
			self._nfft=vel.nvel()/pos.natoms()
			self._natomucell=eig.natomucell()
			self._nrep=self._nfft*pos.natoms()/self._natomucell

	def autoCorr(x):
		result = np.correlate(x, x, mode='full')
		return result[result.size/2:]

	def spctEnrg(self,kindex,mode):

		spatial=2*np.complex(0,1)*np.pi*(self._pos.x()/self._pos.lx()*self._kpt[kindex,0]+
			self._pos.y()/self._pos.ly()*self._kpt[kindex,1]+
			self._pos.z()/self._pos.lz()*self._kpt[kindex,2])
		qdot=np.exp(spatial)*(np.tile(self._eig.x(kindex,mode),(self._nrep,1)).flatten(1)*self._vel.vxts()+
		np.tile(self._eig.y(kindex,mode),(self._nrep,1)).flatten(1)*self._vel.vyts()+
		np.tile(self._eig.z(kindex,mode),(self._nrep,1)).flatten(1)*self._vel.vzts())
		return np.fft.fft(np.correlate(qdot,qdot,mode='full'))

