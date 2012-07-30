## ## ## Kevin Parrish - 7/26/2012	## ## ##
## ## ## Post Processing for GK			## ## ##

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
from scipy import stats
from scipy.optimize import leastsq
import numpy as np
import operator
import os

## ## ## Standard Constants
con_kb = 1.3806E-23
con_hbar = 1.054E-34
con_s2ps = 1E-12	#seconds to picoseconds

## ## ## LJ Parameters
LJ_eps = 1.67E-21
LJ_sigma = 3.4E-10
LJ_mass = 6.6326E-26
LJ_tau = sp.sqrt((LJ_mass*(LJ_sigma**2))/LJ_eps)

## ## ## Green Kubo Simulation Parameters
GK_scaleJ = (LJ_eps)/((LJ_sigma**2)*LJ_tau)	##includes factor of 1/V
GK_s = 5
GK_p = 25000
GK_d = GK_p * GK_s
GK_dt_LJ = 0.002
GK_dt = GK_dt_LJ * LJ_tau
GK_sample_rate = 5
GK_total_steps = 1000000

GK_seeds = [22222, 33333, 44444, 55555, 66666, 77777, 88888, 99999, 10101010, 11111111]
GK_temps = [20, 50, 80, 100]
GK_strains = [-.12, -.10, -.08, -.06, -.04, -.02, -.01, 0, .005, .01, .015, .02, .04, .06, .08, .10, .12]
##### END GK params

## ## ## Mechanics Simulation Parameters
Me_p = 125000
Me_total_steps = GK_total_steps
Me_d = Me_total_steps / Me_p
##### END Me params

GK_num_seeds = len(GK_seeds)
GK_num_temps = len(GK_temps)
GK_num_strains = len(GK_strains)

## ## ## Load stress from file
Me_stress = np.load('post.stress.npy')
			
## ## ## Load kappa from file
kappa = np.load('post.kappa.npy')

## ## ## Load young's modulus from file
young = np.load('post.young.npy')

## ## ## Load phase from file
phase = np.loadtxt('post.phase.txt', dtype=str)

## ## ## Load volume from file
GK_volume = np.load('post.volume.npy')


## ## ## Plot kappa vs strain
artists = []
labels = []
for istrain in range(GK_num_strains):
	for itemp, i in zip(range(GK_num_temps), ['r','g','c','m']):
		if phase[itemp][istrain] == 's' or phase[itemp][istrain] == 'l':
			if phase[itemp][istrain] == 's':
				mk = 'o'
			else:
				mk = 'd'
			plt.semilogy(GK_strains[istrain], kappa[itemp,istrain], c=i, ls='None', marker=mk)
			#label section
			tlabel = 'T = '+ str(GK_temps[itemp])+ 'K, phase = '+ phase[itemp][istrain]
			if tlabel not in labels:
				artists.append(plt.scatter([-1,-1], [0, 0], c=i, marker=mk))
				labels.append(tlabel)

h1 = sorted(zip(artists, labels), key=operator.itemgetter(1)) ## sort handles and labels
artists2, labels2 = zip(*h1)
ll = plt.legend(artists2, labels2, loc='upper right')
for l in ll.get_texts(): ## reduce legend txt size
	l.set_fontsize('x-small')
lt = plt.title('Thermal Conductivities vs Strain')
ly = plt.ylabel('$\kappa (W/m-K)$')
lx = plt.xlabel('$\epsilon$')	
lxm = plt.xlim(-0.15, 0.15)
plt.savefig(('pic.kappa.strain.semi.png'))


## ## ## Fit kappa vs strain
loweps = 0
higheps = 8

x = np.array([-.12,-.10,-.08,-.06,-.04,-.02,-.01,0])

def fitfunc(p, x):
    return (p[0]*np.exp(-x*p[1]) + p[2])

def errfunc(p, x, y):
    return fitfunc(p,x) - y

kparam = np.zeros( (GK_num_temps, 3), dtype=float)
p0 = np.array([1, 20, 1])
for itemp in range(GK_num_temps):
	y = kappa[itemp, loweps:higheps]
	kparam[itemp, :], C, info, msg, success = leastsq(errfunc, p0, args=(x, y), full_output=1)
for itemp, i in zip(range(GK_num_temps), ['r','g','c','m']):
	plt.plot(x, fitfunc(kparam[itemp,:], x), c=i)
	print kparam[itemp, :]
plt.savefig(('pic.kappa.semi.fit.png'))
np.savetxt('post.kparam.txt',kparam)
plt.clf()


## ## ## Plot kappa vs stress semilogy
artists = []
labels = []
for istrain in range(GK_num_strains):
	for itemp, i in zip(range(GK_num_temps), ['r','g','c','m']):
		if phase[itemp][istrain] == 's' or phase[itemp][istrain] == 'l':
			if phase[itemp][istrain] == 's':
				mk = 'o'
			else:
				mk = 'd'
			plt.semilogy(Me_stress[itemp, istrain], kappa[itemp,istrain], c=i, ls='None', marker=mk)
			#label section
			tlabel = 'T = '+ str(GK_temps[itemp])+ 'K, phase = '+ phase[itemp][istrain]
			if tlabel not in labels:
				artists.append(plt.scatter([-1,-1], [0, 0], c=i, marker=mk))
				labels.append(tlabel)

h1 = sorted(zip(artists, labels), key=operator.itemgetter(1)) ## sort handles and labels
artists2, labels2 = zip(*h1)
ll = plt.legend(artists2, labels2, loc='lower left')
for l in ll.get_texts(): ## reduce legend txt size
	l.set_fontsize('x-small')
lt = plt.title('Thermal Conductivities vs Stress')
ly = plt.ylabel('$\kappa (W/m-K)$')
lx = plt.xlabel('$\sigma (Pa)$')	##add units
#lxm = plt.xlim(xmax=4e8)
plt.savefig(('pic.kappa.stress.semi.png'))
plt.clf()


## ## ## Find bulk modulus
def hook(B, x):
	return B * x

def residual(B, x, P):
	return P - hook(B, x)

bulk = np.zeros( (GK_num_temps) )
loweps = 7	#index of low strain inclusive
higheps = 10	#index of high strain exclusive
B0 = 100000	#initial guess

# Create initialV and deltaV for use
initialV = np.zeros( (GK_num_temps), dtype=float)
deltaV = np.zeros( (GK_num_temps, GK_num_strains), dtype=float)
for itemp in range(GK_num_temps):
	initialV[itemp] = GK_volume[itemp, loweps]
	deltaV[itemp,:] = GK_volume[itemp,:] - initialV[itemp]

# Fit bulk modulus
for itemp in range(GK_num_temps):
	x = deltaV[itemp,loweps:higheps] / initialV[itemp]
	P = Me_stress[itemp, loweps:higheps]
	bulk[itemp],cov,infodict,mesg,ier = leastsq(residual, B0, args=(x, P), full_output=True)


## ## ## Plot pressure/volume curve const temp
artists = []
labels = []
for istrain in range(GK_num_strains):
	for itemp, i in zip(range(GK_num_temps), ['r','g','c','m']):
		if phase[itemp][istrain] == 's':
			mk = 'o'
		elif phase[itemp][istrain] == 'l':
			mk = 'd'
		else:
			mk = '2'

		plt.scatter( (deltaV[itemp, istrain] / initialV[itemp]), Me_stress[itemp, istrain], c=i, marker=mk )
		#Label section
		tlabel = 'T = '+ str(GK_temps[itemp])+ 'K, phase = '+ phase[itemp][istrain]
		if tlabel not in labels:
			artists.append(plt.scatter([-4,-4], [4, 4], c=i, marker=mk))
			labels.append(tlabel)

h1 = sorted(zip(artists, labels), key=operator.itemgetter(1)) ## sort handles and labels
artists2, labels2 = zip(*h1)
ll = plt.legend(artists2, labels2, loc='lower right')
for l in ll.get_texts(): ## reduce legend txt size
	l.set_fontsize('x-small')
lt = plt.title('Stress/Strain Curve')
ly = plt.ylabel('$\sigma (Pa)$') 
#ly = plt.ylabel('$\sigma$')
lx = plt.xlabel(r'$\frac{V_c}{V_i}$')	##add units
lxm = plt.xlim(-0.25, 0.25)
plt.savefig(('pic.bulk.1.png'))
lxm = plt.xlim(-0.10, 0.10)
lym = plt.ylim(-1.2e9, 0.3e9)
plt.savefig(('pic.bulk.2.png'))

for itemp, i in zip(range(young[:].size), ['r','g','c','m']):
	plt.plot((deltaV[itemp,loweps:higheps]/initialV[itemp]), hook(bulk[itemp], (deltaV[itemp,loweps:higheps] / initialV[itemp])), c=i, label='For T = '+ str(GK_temps[itemp])+ ', B = %.3g' % bulk[itemp])
ll = plt.legend(loc='upper left')
for l in ll.get_texts(): ## reduce legend txt size
	l.set_fontsize('x-small')
lxm = plt.xlim(-0.005, 0.04)
lym = plt.ylim(-5e7, 2e8)
plt.savefig(('pic.bulk.3.png'))

## ## ## Plot stress/strain curve const temp
artists = []
labels = []
for istrain in range(GK_num_strains):
	for itemp, i in zip(range(GK_num_temps), ['r','g','c','m']):
		if phase[itemp][istrain] == 's':
			mk = 'o'
		elif phase[itemp][istrain] == 'l':
			mk = 'd'
		else:
			mk = '2'
		plt.scatter(GK_strains[istrain], Me_stress[itemp, istrain], c=i, marker=mk)
		#label section
		tlabel = 'T = '+ str(GK_temps[itemp])+ 'K, phase = '+ phase[itemp][istrain]
		if tlabel not in labels:
			artists.append(plt.scatter([-4,-4], [4, 4], c=i, marker=mk))
			labels.append(tlabel)

h1 = sorted(zip(artists, labels), key=operator.itemgetter(1)) ## sort handles and labels
artists2, labels2 = zip(*h1)
ll = plt.legend(artists2, labels2, loc='lower right')
for l in ll.get_texts(): ## reduce legend txt size
	l.set_fontsize('x-small')
lt = plt.title('Stress/Strain Curve')
ly = plt.ylabel('$\sigma (Pa)$')
lx = plt.xlabel('$\epsilon$')	##add units
ll = plt.legend(loc='lower right')
lxm = plt.xlim(-0.15, 0.15)
plt.savefig(('pic.stress.1.png'))
lxm = plt.xlim(-0.09, 0.15)
lym = plt.ylim(-2e9, 0.3e9)
plt.savefig(('pic.stress.2.png'))
lxm = plt.xlim(-0.05, 0.05)
lym = plt.ylim(-6e8, 2e8)
plt.savefig(('pic.stress.3.png'))


"""
## ## ## Overlay with young's modulus fits
def hook(m, x):
	return m*x

higheps = 10
loweps = 7
x = np.zeros( (higheps - loweps), dtype=float)
x = GK_strains[loweps:higheps]
y = np.zeros((5))
x2 = np.zeros((5))
for itemp, i in zip(range(young[:].size), ['r','g','c','m']):
	for j in range(higheps - loweps):
		y[j] = hook(young[itemp], x[j])
		x2[j] = x[j]
	plt.plot(x2, y, c=i, label='Young\'s modulus is %3.2g at temperature of ' % young[itemp] + str(GK_temps[itemp])+ 'K')
ll = plt.legend(loc='upper left')
for l in ll.get_texts(): ## reduce legend txt size
	l.set_fontsize('x-small')
lxm = plt.xlim(-0.0005, 0.0225)
lym = plt.ylim(-0.5e8, 2.5e8)
plt.savefig(('pic.stress.4.png'))
plt.clf()
""" and False





