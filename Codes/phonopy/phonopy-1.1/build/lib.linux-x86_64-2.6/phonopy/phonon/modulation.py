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
import sys
from phonopy.structure.atoms import Atoms
from phonopy.interface.vasp import write_vasp


def write_modulations( dynamical_matrix,
                       primitive,
                       setting,
                       filename="MPOSCAR" ):
    """
    write_modulations:
      ``setting`` is given as sets of band index and amplitude.
      For example:
      setting['q'] = [ 0, 0, 0 ]
      setting['dimension'] = [ 1, 1, 1 ]
      setting['modulations'] = [ [ 1, 2 ], [2, 2.5 ] ]
    """
    
    dimension = setting['dimension']
    argument = setting['argument']
    q = setting['q']
    dm = dynamical_matrix
    dm.set_dynamical_matrix( q )
    eig, eigvec = np.linalg.eigh( dm.get_dynamical_matrix() )
    modulations = []

    for i, ( band_index, amplitude ) in enumerate( setting['modulations'] ):
        modulation = get_modulation( band_index - 1,
                                     eigvec,
                                     primitive,
                                     q,
                                     dimension ) * amplitude 
        cell = get_cell_with_modulation( modulation, primitive, dimension, argument )
        write_vasp( ( filename+"-%03d" ) % (i+1), cell, direct=True )
        
        modulations.append( modulation )

    sum_of_modulations = np.sum( modulations, axis=0 )
    cell = get_cell_with_modulation( sum_of_modulations, primitive, dimension, argument )
    write_vasp( filename, cell, direct=True )

    no_modulations = np.zeros( sum_of_modulations.shape, dtype=complex )
    cell = get_cell_with_modulation( no_modulations, primitive, dimension, argument )
    write_vasp( filename+"-orig", cell, direct=True )

def get_cell_with_modulation( modulation, primitive, dimension, argument ):
    scaled_positions = []
    masses = []
    symbols = []
    for a in range( dimension[0] ):
        for b in range( dimension[1] ):
            for c in range( dimension[2] ):
                for i in range( primitive.get_number_of_atoms() ):
                    p = primitive.get_scaled_positions()[i]
                    scaled_positions.append( p + np.array( [ a,b,c ] ) )
                    masses.append( primitive.get_masses()[i] )
                    symbols.append( primitive.get_chemical_symbols()[i] )

    phase_shift = np.exp( 1j * argument / 180 * np.pi )
    lattice = np.dot( np.diag( dimension ), primitive.get_cell() )
    positions = np.dot( scaled_positions, primitive.get_cell() ) + \
        ( modulation * phase_shift ).imag

    scaled_positions_new = np.dot( positions, np.linalg.inv( lattice ) )
    for p in scaled_positions_new:
        p -= np.floor(p)

    cell = Atoms( cell=lattice,
                  scaled_positions=scaled_positions_new,
                  masses=masses,
                  symbols=symbols,
                  pbc=True )

    return cell
        
def get_modulation( band_index, eigvec, primitive, q, dimension ):
    m = primitive.get_masses()
    r = primitive.get_scaled_positions()
    u = []
    for a in range( dimension[0] ):
        for b in range( dimension[1] ):
            for c in range( dimension[2] ):
                for i, e in enumerate( eigvec[ :, band_index ] ):
                    phase = 2j * np.pi * np.dot( r[ i//3 ] + np.array([a,b,c]), q )
                    u.append( e / np.sqrt( m[ i//3 ] ) * np.exp( phase ) ) 

    return np.array( u ).reshape( -1, 3 )
