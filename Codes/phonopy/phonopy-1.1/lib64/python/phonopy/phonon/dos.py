# Copyright (C) 2011 Atsushi Togo
#
# This file is part of phonopy.
#
# Phonopy is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Phonopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with phonopy.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from phonopy.units import VaspToTHz

class NormalDistribution:
    def __init__(self, sigma):
        self.sigma = sigma

    def calc(self, x):
        return 1.0 / np.sqrt(2 * np.pi) / self.sigma * \
            np.exp(-x**2 / 2.0 / self.sigma**2)
    

class CauchyDistribution:
    def __init__(self, gamma):
        self.gamma = gamma

    def calc(self, x):
        return self.gamma / np.pi / ( x**2 + self.gamma**2 )


class Dos:
    def __init__(self, eigenvalues, weights, factor=VaspToTHz, sigma=None):
        self.eigenvalues = eigenvalues
        self.weights = weights
        if sigma == None:
            self.sigma = factor/100
        else:
            self.sigma = sigma
        self.factor = factor
        self.__frequencies()
        self.set_draw_area()
        # Default smearing 
        self.set_smearing_function('Normal')

    def set_smearing_function(self, function_name):
        """
        function_name ==
        'Normal': smearing is done by normal distribution.
        'Cauchy': smearing is done by Cauchy distribution.
        """
        if function_name == 'Cauchy':
            self.smearing_function = CauchyDistribution(self.sigma)
        else:
            self.smearing_function = NormalDistribution(self.sigma)

    def set_sigma(self, sigma):
        self.sigma = sigma

    def __frequencies(self):
        frequencies = []
        for eig in self.eigenvalues:
            freq = []
            for elem in eig:
                if elem < 0:
                    freq.append(-np.sqrt(-elem))
                else:
                    freq.append(np.sqrt(elem))
            frequencies.append(np.array(freq))

        self.frequencies = np.array(frequencies) * self.factor

    def set_draw_area(self, omega_min=None,
                      omega_max=None, omega_pitch=None):

        if omega_pitch == None:
            self.omega_pitch = \
                float(self.frequencies.max() - self.frequencies.min()) / 200
        else:
            self.omega_pitch = omega_pitch

        if omega_min == None:
            self.omega_min = self.frequencies.min() - self.sigma * 10
        else:
            self.omega_min = omega_min

        if omega_max == None:
            self.omega_max = self.frequencies.max() + self.sigma * 10
        else:
            self.omega_max = omega_max
                    

class TotalDos(Dos):
    def __init__(self, eigenvalues, weights, factor=VaspToTHz, sigma=None):
        Dos.__init__(self, eigenvalues, weights, factor, sigma)

    def get_density_of_states_at_omega(self, omega):
        return np.sum(np.dot(self.weights,
                self.smearing_function.calc(self.frequencies - omega))
                      ) /  np.sum(self.weights)

    def plot_dos(self):
        import matplotlib.pyplot as plt
        plt.plot(self.omegas, self.dos, 'r-')
        plt.grid(True)
        plt.xlabel('Frequency')
        plt.ylabel('Density of states')
        
        return plt
    
    def calculate(self):
        omega = self.omega_min
        dos = []
        while omega < self.omega_max + self.omega_pitch/10 :
            dos.append([omega, self.get_density_of_states_at_omega(omega)])
            omega += self.omega_pitch

        dos = np.array(dos)
        self.omegas = dos[:,0]
        self.dos = dos[:,1]

    def get_dos( self ):
        """
        Return omegas and total dos
        """
        return self.omegas, self.dos

    def write(self):
        file = open('total_dos.dat', 'w')
        file.write("# Sigma = %f\n" % self.sigma)
        for omega, dos in zip( self.omegas, self.dos ):
            file.write("%20.10f%20.10f" % (omega, dos))
            file.write("\n")

class PartialDos(Dos):
    def __init__(self, eigenvalues, weights, eigenvectors,
                 factor=VaspToTHz, sigma=None):
        Dos.__init__(self, eigenvalues, weights, factor, sigma)
        self.eigenvectors2 = (np.abs(eigenvectors))**2

    def get_partial_density_of_states_at_omega(self, index, amplitudes):
        eigvec2 = self.eigenvectors2
        return np.sum(np.dot(self.weights, eigvec2[:,index,:] * amplitudes )
                      ) / np.sum(self.weights)

    def calculate(self):
        omega = self.omega_min
        pdos = []
        omegas = []
        while omega < self.omega_max + self.omega_pitch/10 :
            omegas.append(omega)
            axis_dos = []
            amplitudes = self.smearing_function.calc(self.frequencies - omega)
            for i in range(self.eigenvalues.shape[1]):
                axis_dos.append(self.get_partial_density_of_states_at_omega(i, amplitudes))
            omega += self.omega_pitch
            pdos.append(axis_dos)

        self.partial_dos = np.array(pdos).transpose()
        self.omegas = np.array(omegas)

    def get_partial_dos(self):
        """
        omegas: Sampling frequencies
        partial_dos:
          [[elem1-freq1, elem1-freq2, ... ],
           [elem2-freq1, elem2-freq2, ... ],
           ... ]

          where
           elem1: atom1-x compornent
           elem2: atom1-y compornent
           elem3: atom1-z compornent
           elem4: atom2-x compornent
           ...
        """
        return self.omegas, self.partial_dos

    def plot_pdos(self, indices):
        import matplotlib.pyplot as plt
        plt.grid(True)
        plt.xlim( self.omega_min, self.omega_max )
        plt.xlabel('Frequency')
        plt.ylabel('Partial density of states')
        plots = []

        num_atom = self.eigenvalues.shape[1]/3

        if indices == None:
            indices = []
            for i in range(num_atom):
                indices.append([i+1])

        for set_for_sum in indices:
            pdos_sum = np.zeros(self.omegas.shape, dtype=float)
            for i in set_for_sum:
                if i > num_atom or i < 1:
                    print "Your specified atom number is out of range."
                    raise ValueError
                pdos_sum += self.partial_dos[(i-1)*3:i*3].sum(axis=0)
            plots.append(plt.plot(self.omegas, pdos_sum))

        plt.legend(plots, indices)

        return plt


    def write(self):
        file = open('partial_dos.dat', 'w')
        file.write("# Sigma = %f\n" % self.sigma)
        num_mode = self.eigenvalues.shape[1]
        for omega, pdos in zip(self.omegas, self.partial_dos.transpose()):
            file.write("%20.10f" % omega)
            file.write(("%20.10f" * num_mode) % tuple(pdos))
            file.write("\n")



