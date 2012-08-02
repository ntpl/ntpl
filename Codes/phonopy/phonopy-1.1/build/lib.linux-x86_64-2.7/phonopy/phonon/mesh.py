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
from phonopy.structure.symmetry import get_ir_reciprocal_mesh
from phonopy.units import VaspToTHz

class Mesh:
    def __init__( self,
                  dynamical_matrix,
                  cell,
                  mesh,
                  shift=[0.,0.,0.],
                  is_time_reversal=False,
                  is_symmetry=True,
                  is_eigenvectors=False,
                  factor=VaspToTHz,
                  symprec=1e-5 ):
        self.mesh = np.array(mesh)
        self.shift = np.array(shift)
        self.is_time_reversal = is_time_reversal
        self.is_symmetry = is_symmetry
        self.is_eigenvectors = is_eigenvectors
        self.factor = factor
        self.cell = cell
        self.symprec = symprec
        self.dynamical_matrix = dynamical_matrix

        self.qpoints = None
        self.weights = None
        self.__qpoints()

        self.eigenvalues = None
        self.eigenvectors = None
        self.__eigenvalues()

    def __qpoints(self):
        shift = np.array(self.shift)
        diff = np.abs(shift * 2 - (shift * 2).round())
        if ( diff < self.symprec ).all() and self.is_symmetry:
            diff = np.abs(shift - shift.round())
            self.__qpoints_symmetry(np.logical_xor(
                    (diff > self.symprec),
                    (self.mesh % 2 == 0)))
        else:
            self.__qpoints_no_symmetry()

    def __qpoints_symmetry(self, is_shift):
        map, grid = get_ir_reciprocal_mesh(self.mesh, self.cell,
                                           is_shift * 1,
                                           self.is_time_reversal,
                                           self.symprec)
        ir_list = np.unique(map)
        weights = np.zeros(ir_list.shape[0], dtype=int)
        qpoints = np.zeros((ir_list.shape[0], 3), dtype=float)

        for i, g in enumerate(ir_list):
            weights[i] = np.sum(map == g)
            qpoints[i] = ( grid[g] + is_shift * 0.5 ) / self.mesh
            qpoints[i] -= (qpoints[i] > 0.5) * 1

        self.qpoints = qpoints
        self.weights = weights

    def get_unique_points(self, map):
        ir_list = []
        for x in map:
            is_found = True
            for ir in ir_list:
                if x == ir:
                    is_found = False
                    break

            if is_found:
                ir_list.append(x)

        return np.array(ir_list)

    def __qpoints_no_symmetry(self):
        qpoints = []
        mesh = self.mesh
        shift = self.shift / mesh
        for i in (0, 1, 2):
            if mesh[i] % 2 == 0:
                shift[i] += 0.5 / mesh[i]
                
        for i in range(mesh[2]):
            for j in range(mesh[1]):
                for k in range(mesh[0]):
                    q = np.array([ float(k) / mesh[0],
                                   float(j) / mesh[1],
                                   float(i) / mesh[2] ]) + shift
                    qpoints.append(np.array([ q[0] - ( q[0] > 0.5 ),
                                              q[1] - ( q[1] > 0.5 ),
                                              q[2] - ( q[2] > 0.5 ) ]))

        self.qpoints = np.array(qpoints)
        self.weights = np.ones(self.qpoints.shape[0], dtype=int)

    def get_qpoints(self):
        return self.qpoints

    def get_weights(self):
        return self.weights

    def __eigenvalues( self ):
        eigs = []
        vecs = []
        for q in self.qpoints:
            self.dynamical_matrix.set_dynamical_matrix( q )
            dm = self.dynamical_matrix.get_dynamical_matrix()

            if self.is_eigenvectors:
                val, vec = np.linalg.eigh( dm )
                eigs.append( val.real )
                vecs.append( vec )
            else:
                eigs.append( np.linalg.eigvalsh( dm ).real )

        self.eigenvalues = np.array(eigs)
        if self.is_eigenvectors:
            self.eigenvectors = np.array(vecs)

    def get_eigenvalues(self):
        return self.eigenvalues

    def get_eigenvectors(self):
        """
        Eigenvectors is a numpy array of three dimension.
        The first index runs through q-points.
        In the second and third indices, eigenvectors obtained
        using numpy.linalg.eigh are stored.
        
        The third index corresponds to the eigenvalue's index.
        The second index is for atoms [ x1, y1, z1, x2, y2, z2, ... ].
        """
        return self.eigenvectors


    def write_yaml(self):
        file = open('mesh.yaml', 'w')
        eigenvalues = self.eigenvalues
        natom = self.cell.get_number_of_atoms()
        file.write("mesh: [ %5d, %5d, %5d ]\n" % tuple( self.mesh ) )
        file.write("nqpoint: %-7d\n" % self.qpoints.shape[0])
        file.write("natom:   %-7d\n" % natom )
        file.write("phonon:\n")

        for i, q in enumerate(self.qpoints):
            file.write("- q-position: [ %12.7f, %12.7f, %12.7f ]\n" % tuple(q))
            file.write("  weight: %-5d\n" % self.weights[i])
            file.write("  band:\n")

            for j, eig in enumerate(eigenvalues[i]):
                file.write("  - # %d\n" % ( j+1 ))
                if eig < 0:
                    freq = -np.sqrt(-eig)
                else:
                    freq = np.sqrt(eig)
                file.write("    frequency:  %15.10f\n" % (freq * self.factor))

                if self.is_eigenvectors:
                    file.write("    eigenvector:\n")
                    for k in range(natom):
                        file.write("    - # atom %d\n" % (k+1))
                        for l in (0,1,2):
                            file.write("      - [ %17.14f, %17.14f ]\n" %
                                       (self.eigenvectors[i,k*3+l,j].real,
                                        self.eigenvectors[i,k*3+l,j].imag))
            file.write("\n")
