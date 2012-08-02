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
from phonopy.structure.cells import get_reduced_bases
from phonopy.harmonic.dynamical_matrix import get_equivalent_smallest_vectors

def get_force_constants( set_of_forces,
                         symmetry,
                         supercell,
                         is_tensor_symmetry=False ):

    force_constants = run_force_constants( supercell,
                                           symmetry,
                                           set_of_forces )

    if is_tensor_symmetry:
        set_tensor_symmetry( force_constants, supercell, symmetry )

    return force_constants

def symmetrize_force_constants( force_constants, iteration=3 ):
    for i in range( iteration ):
        set_permutation_symmetry( force_constants )
        set_translational_invariance( force_constants )

def run_force_constants( supercell, symmetry, set_of_forces ):
    """
    Bare force_constants is returned.

    Force constants, Phi, are calculated from sets for forces, F, and
    atomic displacement, d:
      Phi = -F / d
    This is solved by matrix pseudo-inversion.
    Crsytal symmetry is included when creating F and d matrices.

    force_constants[ i, j, a, b ]
      i: Atom index of finitely displaced atom.
      j: Atom index at which force on the atom is measured.
      a, b: Cartesian direction indices = (0, 1, 2) for i and j, respectively
    """
    
    force_constants = np.zeros( ( supercell.get_number_of_atoms(),
                                  supercell.get_number_of_atoms(),
                                  3, 3), dtype=float)

    # Fill force_constants[ displaced_atoms, all_atoms_in_supercell ]
    atom_list_done = get_force_constant_disps( force_constants,
                                               supercell,
                                               set_of_forces,
                                               symmetry )

    # Distribute non-equivalent force constants to those equivalent
    for i in range( supercell.get_number_of_atoms() ):
        if not ( i in atom_list_done ):
            distribute_force_constant( force_constants,
                                       i,
                                       atom_list_done,
                                       supercell,
                                       symmetry )

    return force_constants

def distribute_force_constant( force_constants,
                               atom_number,
                               atom_list_done,
                               supercell,
                               symmetry,
                               is_symmetrize=False ):

    positions = supercell.get_scaled_positions()
    symprec = symmetry.get_symmetry_tolerance()
    rotations = symmetry.get_symmetry_operations()['rotations']
    trans = symmetry.get_symmetry_operations()['translations']

    map_atom_disps, map_syms = \
        get_atomic_mapping_by_symmetry( atom_list_done,
                                        atom_number,
                                        rotations, 
                                        trans,
                                        positions,
                                        symprec )

    if not is_symmetrize:
        map_atom_disps = [ map_atom_disps[0], ]
        map_syms = [ map_syms[0], ]

    for i, pos_i in enumerate( positions ):
        for map_atom_disp, map_sym in zip( map_atom_disps, map_syms ):
            rot_pos = np.dot( pos_i, rotations[ map_sym ].T ) + trans[ map_sym ]
            rot_atom = -1
            for j, pos_j in enumerate( positions ):
                diff = pos_j - rot_pos
                if ( abs( diff - diff.round() ) < symprec ).all():
                    rot_atom = j
                    break
    
            if rot_atom < 0:
                print 'Input forces are not enough to calculate force constants,'
                print 'or something wrong (e.g. crystal structure does not match).'
                raise ValueError
    
            # L R L^-1
            rot_cartesian = similarity_transformation( supercell.get_cell().T,
                                                       rotations[ map_sym ] )
    
            # R^-1 P R (inverse transformation)
            force_constants[ atom_number, i ] += \
                similarity_transformation( rot_cartesian.T,
                                           force_constants[ map_atom_disp, rot_atom ] )

        force_constants[ atom_number, i ] /= len( map_atom_disps )
    
def get_atomic_mapping_by_symmetry( atom_list_done,
                                    atom_number,
                                    rotations,
                                    translations,
                                    positions,
                                    symprec=1e-5 ):
    """
    The mapping from an atom to the atom in the atom list done.
    """

    map_atoms = []
    map_syms = []
    for i, (r, t) in enumerate( zip( rotations, translations ) ):
        rot_pos = np.dot( positions[atom_number], r.T ) + t
        for j in atom_list_done:
            diff = positions[j] - rot_pos
            if ( abs( diff -diff.round() ) < symprec ).all():
                map_atoms.append( j )
                map_syms.append( i )
                break

    if len( map_atoms ) == 0:
        print 'Input forces are not enough to calculate force constants,'
        print 'or something wrong (e.g. crystal structure does not match).'
        raise ValueError

    return map_atoms, map_syms

def get_force_constant_disps( force_constants, supercell, set_of_forces, symmetry ):
    """
    Phi = -F / d
    """
    
    disp_atom_list = np.unique([ forces.get_atom_number() for forces in set_of_forces ])

    for disp_atom_number in disp_atom_list:
        site_symmetry = symmetry.get_site_symmetry(disp_atom_number)
        symmetry_matrices = get_symmetry_matrices(supercell, site_symmetry)
        rot_disps = []
        row_forces = []

        # Bind several displacements of a displaced atom
        # with symmetry operations
        for forces in set_of_forces:
            if forces.get_atom_number() != disp_atom_number:
                continue

            displacement = forces.get_displacement()
            # Displacement * Rotation (U * A)
            rot_disps.append( get_rotated_displacements( displacement,
                                                         symmetry_matrices ) )
            row_forces.append( forces.get_forces() )

        # Bind forces for several displacements and symmetry
        # operations of a displaced atom
        rot_disps = np.array(rot_disps).reshape(-1, 9)
        inv = np.linalg.pinv(rot_disps)
        for i in range(supercell.get_number_of_atoms()):
            combined_forces = []
            for forces in row_forces:
                # Combined forces (F)
                combined_forces.append(get_combined_force(supercell,
                                                          i,
                                                          disp_atom_number,
                                                          forces,
                                                          site_symmetry,
                                                          symmetry.get_symmetry_tolerance()))

            combined_forces = np.array(combined_forces).reshape(-1, 1)
            # Force constant (P) = -(U * A)^-1 * (F)
            force_constants[disp_atom_number,i] = \
                -np.dot( inv, combined_forces ).reshape(3, 3)

    return disp_atom_list

def set_tensor_symmetry( force_constants, supercell, symmetry ):
    """
    Full force constants are symmetrized using crystal symmetry.
    This method extracts symmetrically equivalent sets of atomic pairs and
    take sum of their force constants and average the sum.
    
    Since get_force_constant_disps may include crystal symmetry, this method
    is usually meaningless.
    """

    positions = supercell.get_scaled_positions()
    symprec = symmetry.get_symmetry_tolerance()
    rotations = symmetry.get_symmetry_operations()['rotations']
    translations = symmetry.get_symmetry_operations()['translations']

    fc_bak = force_constants.copy()

    # Create mapping table between an atom and the symmetry operated atom
    # map[ i, j ]
    # i: atom index
    # j: operation index
    mapping = []
    for pos_i in positions:
        map_local = []
        for rot, trans in zip( rotations, translations ):
            rot_pos = np.dot( pos_i, rot.transpose() ) + trans
            for j, pos_j in enumerate( positions ):
                diff = pos_j - rot_pos
                if ( abs( diff -diff.round() ) < symprec ).all():
                    map_local.append( j )
                    break
        mapping.append( map_local )
    mapping = np.array( mapping )

    # Look for the symmetrically equivalent force constant tensors
    for i, pos_i in enumerate( positions ):
        for j, pos_j in enumerate( positions ):
            tmp_fc = np.zeros( ( 3, 3 ), dtype=float )
            for k, rot in enumerate( rotations ):
                cart_rot = similarity_transformation( supercell.get_cell().transpose(), rot )

                # Reverse rotation of force constant is summed
                tmp_fc += similarity_transformation( cart_rot.transpose(),
                                                     fc_bak[ mapping[ i, k ],
                                                             mapping[ j, k ] ] )
            # Take average and set to new force cosntants
            force_constants[ i, j ] = tmp_fc / len( rotations )

def set_translational_invariance( force_constants ):
    """
    Translational invariance is imposed.  This is quite simple
    implementation, which is just take sum of the force constants in
    an axis and an atom index. The sum has to be zero due to the
    translational invariance. If the sum is not zero, this error is
    uniformly subtracted from force constants.
    """
    for i in range(force_constants.shape[1]):
        for j in range(force_constants.shape[2]):
            for k in range(force_constants.shape[3]):
                force_constants[:,i,j,k] -= \
                    np.sum( force_constants[:,i,j,k] ) / force_constants.shape[0]

def set_permutation_symmetry( force_constants ):
    """
    Phi_ij_ab = Phi_ji_ba
    
    i, j: atom index
    a, b: Cartesian axis index

    This is not necessary for harmonic phonon calculation because this
    condition is imposed when making dynamical matrix Hermite in
    dynamical_matrix.py.
    """
    fc_copy = force_constants.copy()
    for i in range(force_constants.shape[0]):
        for j in range(force_constants.shape[1]):
            force_constants[ i, j ] = ( force_constants[ i, j ] + fc_copy[ j, i ].T ) / 2

def rotational_invariance(force_constants, supercell, primitive, symprec=1e-5):
    """
    *** Under development ***
    Just show how force constant is close to the condition of rotational invariance,
    """
    print "Check rotational invariance ..."

    fc = force_constants 
    p2s = primitive.get_primitive_to_supercell_map()

    abc = "xyz"
    
    for pi, p in enumerate( p2s ):
        for i in range(3):
            mat = np.zeros( ( 3, 3 ), dtype=float )
            for s in range( supercell.get_number_of_atoms() ):
                vecs = np.array( get_equivalent_smallest_vectors(
                        s, p, supercell, primitive.get_cell(), symprec ) )
                m = len( vecs )
                v = np.dot( vecs[:,:].sum(axis=0) / m, primitive.get_cell() )
                for j in range(3):
                    for k in range(3):
                        mat[j,k] += fc[p,s,i,j] * v[k] - fc[p,s,i,k] * v[j]

            print "Atom %d %s" % ( p+1, abc[i] )
            for vec in mat:
                print "%10.5f %10.5f %10.5f" % tuple( vec )

def force_constants_log(force_constants):
    fs = force_constants
    for i, fs_i in enumerate(fs):
        for j, fs_j in enumerate(fs_i):
            for v in fs_j:
                print "force constant (%d - %d): %10.5f %10.5f %10.5f" % (i+1,j+1,v[0],v[1],v[2])


# Shared functions to calculate force constant
def get_symmetry_matrices( cell, site_symmetry ):
    """
    Transformation of 2nd order force constant

    In the phonopy implementation (Cartesian coords.)

    (R.F)^T = -(R.U)^T Psi' --> F^T = -U^T.R^T.Psi'.R
    Therefore,
    Psi = R^T.Psi'.R --> Psi' = R.Psi.R^T

    The symmetrical relation between Psi and Psi' can be represented
    by a 9x9 matrix. What we want is transformation matrix A defined
    by

    P' = A.P

    where P' and P are the 9x1 matrices and A is the 9x9 matrices.
    """
    matrices = []
    for reduced_rot in site_symmetry:
        mat = []
        rot = similarity_transformation( cell.get_cell().T, reduced_rot )
        for i in range( 3 ):
            for j in range( 3 ):
                psi = []
                for k in range( 3 ):
                    for l in range( 3 ):
                        psi.append( rot[ i, k ] * rot[ j, l ] )
                mat.append( psi )
        matrices.append( mat )
    return np.array(matrices)

def similarity_transformation(rot, mat):
    """ R x M x R^-1 """
    return np.dot( rot, np.dot( mat, np.linalg.inv(rot) ) )
                
def get_rotated_displacements(displacement, symmetry_matrices):
    """
    U x A
                                                [123456789]
                                                [2        ]
                                                [3        ]
      [ d_x  0   0   d_y  0   0   d_z  0   0  ] [4        ]
    U [  0  d_x  0    0  d_y  0    0  d_z  0  ] [5   A    ]
      [  0   0  d_x   0   0  d_y   0   0  d_z ] [6        ]
                                                [7        ]
                                                [8        ]
                                                [9        ]
    """
    rot_disps = []
    for sym in symmetry_matrices:
        rot_disps.append( np.dot( expand_displacement(displacement), sym ) )
    return rot_disps

def expand_displacement(displacement):
    """
    [ d_x  0   0   d_y  0   0   d_z  0   0  ]
    [  0  d_x  0    0  d_y  0    0  d_z  0  ]
    [  0   0  d_x   0   0  d_y   0   0  d_z ]
    """
    d = displacement
    disp = np.hstack( ( np.eye(3)*d[0], np.eye(3)*d[1], np.eye(3)*d[2] ) )
    return disp

def get_rotated_atom_number( cell,
                             atom_number,
                             center_atom_number,
                             rotation,
                             symprec=1e-5 ):
    """
    Return atom number sent by site symmetry
    """
    pos = cell.get_scaled_positions()
    atomic_numbers = cell.get_atomic_numbers()
    rot_pos = np.dot( pos[atom_number] - pos[center_atom_number],
                      rotation.T ) + pos[center_atom_number]

    for i, p in enumerate( pos ):
        if not atomic_numbers[i] == atomic_numbers[atom_number]:
            continue

        diff = p - rot_pos
        if ( abs( diff - diff.round() ) < symprec ).all():
            return i

    print "Phonopy enconters symmetry problem."
    return -1
        
def get_combined_force( cell,
                        atom_number,
                        center_atom_number,
                        forces,
                        site_symmetry,
                        symprec=1e-5 ):
    """
    Pack forces on atoms translated by site symmetry
    
    The relation:
    R [ F(r) ] = F( R.r )
    where R is a rotation cetering at displaced atom.
    (This is not the transformation of a function,
     but just the rotation of force vector at r. )
    """
    rot_forces = []
    for sym in site_symmetry:
        rot_atom_number = get_rotated_atom_number( cell,
                                                   atom_number,
                                                   center_atom_number,
                                                   sym,
                                                   symprec)
        rot_forces.append(forces[rot_atom_number])
    return rot_forces

