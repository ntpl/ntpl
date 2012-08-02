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

class BandStructure:
    def __init__( self, bands, dynamical_matrix, cell,
                  is_eigenvectors=False,
                  factor=VaspToTHz,
                  verbose=False ):
        self.dynamical_matrix = dynamical_matrix
        self.cell = cell
        self.factor = factor
        self.is_eigenvectors = is_eigenvectors
        self.bands = bands
        self.distances = []
        self.qpoints = []
        self.distance = 0.
        self.special_point = [0.]
        self.eigenvalues = None
        self.eigenvectors = None
        self.set_band( verbose=verbose )

    def get_distances( self ):
        return self.distances

    def get_qpoints( self ):
        return self.qpoints

    def get_eigenvalues( self ):
        return self.eigenvalues

    def get_eigenvectors( self ):
        return self.eigenvectors
    
    def set_initial_point(self, qpoint):
        self.lastq = qpoint.copy()

    def set_point(self, qpoint):
        self.qpoints.append(qpoint.copy())
        self.distance += np.linalg.norm(
            np.dot( qpoint-self.lastq,
                    np.linalg.inv(self.cell.get_cell()).transpose() )
            )
        self.distances.append(self.distance)
        self.lastq = qpoint.copy()

    def set_band(self, nac=False, verbose=False):
        eigvals = []
        eigvecs = []

        is_nac = self.dynamical_matrix.is_nac()

        for band in self.bands:
            self.set_initial_point(band[0])

            if is_nac:
                # One of end points has to be Gamma point.
                if np.linalg.norm( band[0] ) < 0.0001 or \
                        np.linalg.norm( band[-1] ) < 0.0001:
                    q_direction = band[0] - band[-1]
                else:
                    q_direction = None
            
            for q in band:
                self.set_point(q)
                if is_nac:
                    self.dynamical_matrix.set_dynamical_matrix(
                        q, q_direction, verbose )
                else:
                    self.dynamical_matrix.set_dynamical_matrix(q, verbose)
                dm = self.dynamical_matrix.get_dynamical_matrix()

                if self.is_eigenvectors:
                    val, vec = np.linalg.eigh(dm)
                    eigvals.append(val.real)
                    eigvecs.append(vec)
                else:
                    eigvals.append(np.linalg.eigvalsh(dm).real)

            self.special_point.append(self.distance)

        self.eigenvalues = np.array(eigvals)
        if self.is_eigenvectors:
            self.eigenvectors = np.array(eigvecs)

    def plot_band(self, symbols=None):
        import matplotlib.pyplot as plt

        for eigs in self.eigenvalues.transpose():
            line = []
            for val in eigs:
                if val > 0:
                    line.append(np.sqrt(val))
                else:
                    line.append(-np.sqrt(-val))
            plt.plot(self.distances, np.array(line) * self.factor, 'r-')

        plt.ylabel('Frequency')
        plt.xlabel('Wave vector')
        if (not symbols==None) and len( symbols )==len( self.special_point ):
            plt.xticks(self.special_point, symbols )
        else:
            plt.xticks(self.special_point, [''] * len(self.special_point))
        plt.xlim(0, self.distance)
        plt.axhline(y=0, linestyle=':', linewidth=0.5, color='b')
        return plt

    def write_yaml(self):
        file = open('band.yaml', 'w')
        eigenvalues = self.eigenvalues
        natom = self.cell.get_number_of_atoms()
        file.write("nqpoint: %-7d\n" % len(self.qpoints))
        file.write("npath: %-7d\n" % len(self.bands))
        file.write("natom: %-7d\n" % (natom))
        file.write("phonon:\n")
        for i, q in enumerate(self.qpoints):
            file.write("- q-position: [ %12.7f, %12.7f, %12.7f ]\n" % tuple(q))
            file.write("  distance: %12.7f\n" % self.distances[i])
            file.write("  band:\n")
            for j, eig in enumerate(eigenvalues[i]):
                if eig < 0:
                    freq = -np.sqrt(-eig)
                else:
                    freq = np.sqrt(eig)
                file.write("  - # %d\n" % ( j+1 ))
                file.write("    frequency: %15.10f\n" % (freq * self.factor))

                if self.is_eigenvectors:
                    file.write("    eigenvector:\n")
                    for k in range(natom):
                        file.write("    - # atom %d\n" % (k+1))
                        for l in (0,1,2):
                            file.write("      - [ %17.14f, %17.14f ]\n" %
                                       (self.eigenvectors[i,k*3+l,j].real,
                                        self.eigenvectors[i,k*3+l,j].imag))
            file.write("\n")





