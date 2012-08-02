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
import cmath
from phonopy.units import VaspToTHz

def write_yaml( qpoints,
                cell,
                dynamical_matrix, 
                is_eigenvectors=False,
                factor=VaspToTHz ):
    num_atom = cell.get_number_of_atoms()
    m = cell.get_masses()
    names = cell.get_chemical_symbols()
    positions = cell.get_scaled_positions()
    lattice = cell.get_cell()

    file = open('qpoints.yaml', 'w')
    file.write("nqpoint: %-7d\n" % qpoints.shape[0])
    file.write("natom:   %-7d\n" % num_atom )
    file.write("atom-info:\n")
    for mass, name in zip( m, names ):
        file.write("- { name: %2s, mass: %10.5f }\n" % ( name, mass ))
    
    file.write("real-basis:\n")
    file.write("- [ %20.15f, %20.15f, %20.15f ]\n" % ( tuple( lattice[0] ) ) )
    file.write("- [ %20.15f, %20.15f, %20.15f ]\n" % ( tuple( lattice[1] ) ) )
    file.write("- [ %20.15f, %20.15f, %20.15f ]\n" % ( tuple( lattice[2] ) ) )

    file.write("position:\n")
    for pos in positions:
        file.write("- [ %20.15f, %20.15f, %20.15f ]\n" % ( tuple( pos ) ) )
        

    file.write("phonon:\n")

    for q in qpoints:
        dynamical_matrix.set_dynamical_matrix(q)
        dm = dynamical_matrix.get_dynamical_matrix()

        file.write("- q-position: [ %12.7f, %12.7f, %12.7f ]\n" % tuple(q))
        file.write("  q-point:\n")
            
        if is_eigenvectors:
            eigenvalues, eigenvectors = np.linalg.eigh( dm )
        else:
            eigenvalues = np.linalg.eigvalsh( dm )

        for j, eig in enumerate( eigenvalues ):
            if eig < 0:
                freq = -np.sqrt( -eig )
            else:
                freq = np.sqrt( eig )
            file.write("  - # %d\n" % ( j+1 ) )
            file.write("    frequency: %15.10f\n" % ( freq * factor) )

            if is_eigenvectors:
                file.write("    eigenvector:\n")
                for k in range( num_atom ):
                    file.write("    - # atom %d\n" % ( k+1 ))
                    for l in (0,1,2):
                        file.write("      - [ %17.14f, %17.14f ]\n" %
                                   ( eigenvectors[k*3+l,j].real,
                                     eigenvectors[k*3+l,j].imag ))

                file.write("    eigenvector_time_aligned:\n")
                for k in range( num_atom ):
                    file.write("    - # atom %d, freq*sqrt(m) %f, [%f %f %f]\n" %
                               ( ( k+1, freq * factor * np.sqrt( m[k] ) ) + tuple( positions[k] ) ) )
                    # Phase of dot( q, r ) is added.
                    eig_aligned = eigenvectors[k*3:(k+1)*3, j] * np.exp( 2j * np.pi * np.dot( positions[k], q ) )
                    for l in (0,1,2):
                        file.write("      - [ %17.14f, %17.14f ] # %7.2f, %7.4f\n" %
                                   ( eig_aligned[l].real,
                                     eig_aligned[l].imag,
                                     cmath.phase( eig_aligned[l] ) / np.pi * 180,
                                     abs( eig_aligned[l] ) ) )
        file.write("\n")
