
class Nmd:
	"""
	Creates an nmd object for the nmd calculation.
	nmdpy.Nmd(vel,pos)
	Parameters
	----------
		vel : array
			time series of atomic velocities from MD
		pos : array
			equilibrium atomic positions
	"""
	def __init__(self,
			vel,
			pos):
			self._vel=vel
			self._pos=pos

	def spctEnrg(self,eigvec,k):
		qdot=eigvec.x*self._vel.x+
			eigvec.y*self._vel.y+
			eigvec.y*self._vel.z


