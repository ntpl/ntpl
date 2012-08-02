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
from phonopy.structure.atoms import Atoms

def get_supercell(unitcell, supercell_matrix, symprec = 1e-5):
    """Build supercell from supercell matrix
    In this function, unit cell is considered 
    [1,0,0]
    [0,1,0]
    [0,0,1].
    Supercell matrix is given by relative ratio, e.g,
    [-1, 1, 1]
    [ 1,-1, 1]  is for FCC from simple cubic.
    [ 1, 1,-1]
    In this case multiplicities of surrounding simple lattice are [2,2,2].

    First, create supercell with surrounding simple lattice.
    Second, trim the surrounding supercell with the target lattice.
    """
    
    mat = np.array(supercell_matrix)
    multi = get_surrounding_frame(mat)
    surrounding_cell = get_simple_supercell( multi, unitcell )

    # Trim the simple supercell by the supercell matrix
    trim_frame = np.array( [mat[0] / float(multi[0]),
                            mat[1] / float(multi[1]),
                            mat[2] / float(multi[2])] )
    supercell, mapping = trim_cell( trim_frame,
                                    surrounding_cell,
                                    symprec )
    simple_s2u_map = []
    num_uatom = unitcell.get_number_of_atoms()
    multiplicity = supercell.get_number_of_atoms() / num_uatom
    for i in range( num_uatom ):
        simple_s2u_map += [ i * multiplicity ] * multiplicity
    
    s2u_map = []
    for i in mapping:
        s2u_map.append( simple_s2u_map[ i ] )
    
    return Supercell( supercell, s2u_map )

def get_simple_supercell( multi, unitcell ):
    # Scaled positions within the frame, i.e., create a supercell that
    # is made simply to multiply the input cell.
    positions = unitcell.get_scaled_positions()
    numbers = unitcell.get_atomic_numbers()
    masses = unitcell.get_masses()
    lattice = unitcell.get_cell()
    
    positions_multi = []
    numbers_multi = []
    masses_multi = []
    for l, pos in enumerate(positions):
        for i in range(multi[2]):
            for j in range(multi[1]):
                for k in range(multi[0]):
                    positions_multi.append([ (pos[0] + k) / multi[0],
                                             (pos[1] + j) / multi[1],
                                             (pos[2] + i) / multi[2] ])
                    numbers_multi.append(numbers[l])
                    masses_multi.append(masses[l])

    return Atoms(numbers = numbers_multi,
                 masses = masses_multi,
                 scaled_positions = positions_multi,
                 cell = np.dot( np.diag( multi ), lattice ),
                 pbc=True)

def get_surrounding_frame(supercell_matrix):
    # Build a frame surrounding supercell lattice
    # For example,
    #  [2,0,0]
    #  [0,2,0] is the frame for FCC from simple cubic.
    #  [0,0,2]
    m = np.array( supercell_matrix )
    axes = np.array( [ [ 0, 0, 0],
                       m[:,0],
                       m[:,1],
                       m[:,2],
                       m[:,1] + m[:,2],
                       m[:,2] + m[:,0],
                       m[:,0] + m[:,1],
                       m[:,0] + m[:,1] + m[:,2] ] )
    frame = [ max( axes[:,i] ) - min( axes[:,i] ) for i in (0,1,2) ]
    return frame

def trim_cell(relative_axes, cell, symprec):
    # relative_axes: relative axes to supercell axes
    # Trim positions outside relative axes

    positions = cell.get_scaled_positions()
    numbers = cell.get_atomic_numbers()
    masses = cell.get_masses()
    lattice = cell.get_cell()
    trimed_lattice = np.dot( relative_axes.T, lattice )
    
    trimed_positions = []
    trimed_numbers = []
    trimed_masses = []
    atom_map = []
    for i, pos in enumerate(positions):
        is_overlap = False
        # tmp_pos is the position with respect to the relative axes
        tmp_pos = np.dot(pos, np.linalg.inv(relative_axes).T)
        # Positions are folded into trimed cell and
        # overlap positions are removed.
        for t_pos in trimed_positions:
            diff = t_pos - tmp_pos
            diff -= diff.round()
            if np.linalg.norm( np.dot( diff, trimed_lattice ) ) < symprec:
                is_overlap = True
                break

        if not is_overlap:
            trimed_positions.append(tmp_pos - np.floor(tmp_pos))
            trimed_numbers.append(numbers[i])
            trimed_masses.append(masses[i])
            atom_map.append(i)

    trimed_cell = Atoms( numbers = trimed_numbers,
                         masses = trimed_masses,
                         scaled_positions = trimed_positions,
                         cell = trimed_lattice,
                         pbc=True)

    return trimed_cell, atom_map

def print_cell(cell, mapping=None):
    symbols = cell.get_chemical_symbols()
    masses = cell.get_masses()
    lattice = cell.get_cell()
    print "Lattice vectors:"
    print "  a %20.15f %20.15f %20.15f" % tuple( lattice[0] )
    print "  b %20.15f %20.15f %20.15f" % tuple( lattice[1] )
    print "  c %20.15f %20.15f %20.15f" % tuple( lattice[2] )
    print "Atomic positions (fractional):"
    for i, v in enumerate(cell.get_scaled_positions()):
        print "%5d %-2s%18.14f%18.14f%18.14f %7.3f" % \
            (i+1, symbols[i], v[0], v[1], v[2], masses[i]),
        # print 
        if mapping==None:
            print
        else:
            print ">", mapping[i]+1

def get_supercells_with_displacements( cell, displacements ):
    
    cells = []
    for disp in displacements:
        disp_cartesian = np.dot(disp[1:], cell.get_cell())
        disp_cartesian = disp_cartesian / np.linalg.norm(disp_cartesian) * distance
        positions = cell.get_positions()
        positions[disp[0]] += disp_cartesian
        cells.append( Atoms( numbers = cell.get_atomic_numbers(),
                             masses = cell.get_masses(),
                             positions = positions,
                             cell = cell.get_cell(),
                             pbc = True ) )

    return cells

class Supercell(Atoms):
    def __init__( self, supercell, supercell_to_unitcell_map ):
        Atoms.__init__( self,
                        numbers = supercell.get_atomic_numbers(),
                        masses = supercell.get_masses(),
                        scaled_positions = supercell.get_scaled_positions(),
                        cell = supercell.get_cell(),
                        pbc = True )
        self.s2u_map = supercell_to_unitcell_map

    def get_supercell_to_unitcell_map( self ):
        return self.s2u_map

class Primitive(Atoms):
    def __init__(self, supercell, primitive_frame, symprec = 1e-5):
        """
        primitive_frame ( 3x3 matrix ):
        Primitive lattice is given by
           np.dot( primitive_frame.T, supercell.get_cell())
        """
        self.frame = np.array(primitive_frame)
        self.symprec = symprec
        self.p2s_map = None
        self.s2p_map = None
        self.p2p_map = None
        self.__primitive_cell( supercell )
        self.__supercell_to_primitive_map( supercell.get_scaled_positions() )
        self.__primitive_to_primitive_map()

    def __primitive_cell(self, supercell):
        trimed_cell, p2s_map = trim_cell( self.frame,
                                          supercell,
                                          self.symprec )
        Atoms.__init__( self,
                        numbers = trimed_cell.get_atomic_numbers(),
                        masses = trimed_cell.get_masses(),
                        scaled_positions = trimed_cell.get_scaled_positions(),
                        cell = trimed_cell.get_cell(),
                        pbc = True )

        self.p2s_map = p2s_map

    def get_primitive_to_supercell_map(self):
        return self.p2s_map

    def get_frame(self):
        return self.frame

    def __supercell_to_primitive_map(self, pos):
        inv_F = np.linalg.inv(self.frame)
        s2p_map = []
        for i in range(pos.shape[0]):
            s_pos = np.dot( pos[i], inv_F.T )
            for j in self.p2s_map:
                p_pos = np.dot( pos[j], inv_F.T )
                diff = p_pos - s_pos
                diff -= diff.round()
                if ( abs(diff) < self.symprec ).all():
                    s2p_map.append(j)
                    break
        self.s2p_map = s2p_map

    def get_supercell_to_primitive_map(self):
        return self.s2p_map

    def __primitive_to_primitive_map(self):
        """
        Mapping table from supercell index to primitive index
        in primitive cell
        """
        self.p2p_map = dict([ ( j, i ) for i, j in enumerate( self.p2s_map ) ])

    def get_primitive_to_primitive_map(self):
        return self.p2p_map



#
# Get distance between a pair of atoms
#
def get_distance(atoms, a0, a1, tolerance=1e-5):
    """
    Return the shortest distance between a pair of atoms in PBC
    """
    reduced_bases = get_reduced_bases( atoms.get_cell(), tolerance)
    scaled_pos = np.dot( atoms.get_positions(), np.linalg.inv(reduced_bases) )
    # move scaled atomic positions into -0.5 < r <= 0.5
    for pos in scaled_pos:
        pos -= pos.round()

    # Look for the shortest one in surrounded 3x3x3 cells 
    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append( np.linalg.norm(
                        np.dot(scaled_pos[a0] - scaled_pos[a1] + np.array([i,j,k]),
                               reduced_bases) ) )
    return min(distances)

def get_distance_with_center( atoms, center, atom_num, tolerance=1e-5 ):
    """
    Return the shortest distance to atom from specified position
    """
    reduced_bases = get_reduced_bases( atoms.get_cell(), tolerance)
    scaled_pos = np.dot( atoms.get_positions(), np.linalg.inv(reduced_bases) )
    # move scaled atomic positions into -0.5 < r <= 0.5
    for pos in scaled_pos:
        pos -= pos.round()

    # Look for the shortest one in surrounded 3x3x3 cells 
    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append( np.linalg.norm(
                        np.dot(scaled_pos[atom_num] - center + np.array([i,j,k]),
                               reduced_bases) ) )
    return min(distances)
    

#
# Delaunay reduction
#    
def get_reduced_bases(cell, tolerance=1e-5):
    """
    This is an implementation of Delaunay reduction.
    Some information is found in International table.
    """
    return get_Delaunay_reduction(cell, tolerance)

def get_Delaunay_reduction(lattice, tolerance):
    extended_bases = np.zeros((4,3), dtype=float)
    extended_bases[:3,:] = lattice
    extended_bases[3] = -np.sum(lattice, axis=0)

    for i in range(100):
        if reduce_bases(extended_bases, tolerance):
            break
    if i == 99:
        print "Delaunary reduction is failed."

    shortest = get_shortest_bases_from_extented_bases(extended_bases, tolerance)

    return shortest

def reduce_bases(extended_bases, tolerance):
    metric = np.dot(extended_bases, extended_bases.transpose())
    for i in range(4):
        for j in range(i+1, 4):
            if metric[i][j] > tolerance:
                for k in range(4):
                    if (not k == i) and (not k == j):
                        extended_bases[k] += extended_bases[i]
                extended_bases[i] = -extended_bases[i]
                extended_bases[j] = extended_bases[j]
                return False

    # Reduction is completed.
    # All non diagonal elements of metric tensor is negative.
    return True

def get_shortest_bases_from_extented_bases(extended_bases, tolerance):

    def mycmp(x, y):
        return cmp(np.vdot(x,x), np.vdot(y,y))

    basis = np.zeros((7,3), dtype=float)
    basis[:4] = extended_bases
    basis[4]  = extended_bases[0] + extended_bases[1]
    basis[5]  = extended_bases[1] + extended_bases[2]
    basis[6]  = extended_bases[2] + extended_bases[0]
    # Sort bases by the lengthes (shorter is earlier)
    basis = sorted(basis, cmp=mycmp)
    
    # Choose shortest and linearly independent three bases
    # This algorithm may not be perfect.
    for i in range(7):
        for j in range(i+1, 7):
            for k in range(j+1, 7):
                if abs(np.linalg.det([basis[i],basis[j],basis[k]])) > tolerance:
                    return np.array([basis[i],basis[j],basis[k]])

    print "Delaunary reduction is failed."
    return basis[:3]

#
# Other tiny tools
#    
def get_angles( lattice ):
    a, b, c = get_cell_parameters( lattice )
    alpha = np.arccos(np.vdot(lattice[1], lattice[2]) / b / c ) / np.pi * 180
    beta  = np.arccos(np.vdot(lattice[2], lattice[0]) / c / a ) / np.pi * 180
    gamma = np.arccos(np.vdot(lattice[0], lattice[1]) / a / b ) / np.pi * 180
    return alpha, beta, gamma

def get_cell_parameters( lattice ):
    return np.sqrt( np.dot ( lattice, lattice.transpose() ).diagonal() )

def get_cell_matrix( a, b, c, alpha, beta, gamma ):
    # These follow 'matrix_lattice_init' in matrix.c of GDIS
    alpha *= np.pi / 180
    beta *= np.pi / 180
    gamma *= np.pi / 180
    a1 = a
    a2 = 0.0
    a3 = 0.0
    b1 = np.cos( gamma )
    b2 = np.sin( gamma )
    b3 = 0.0
    c1 = np.cos( beta )
    c2 = ( 2 * np.cos( alpha ) + b1**2 + b2**2 - 2 * b1 * c1 - 1 ) / ( 2 * b2 )
    c3 = np.sqrt( 1 - c1**2 - c2**2 )
    lattice = np.zeros( ( 3, 3 ), dtype=float )
    lattice[ 0, 0 ] = a
    lattice[ 1 ] = np.array( [ b1, b2, b3 ] ) * b
    lattice[ 2 ] = np.array( [ c1, c2, c3 ] ) * c
    return lattice

if __name__ == '__main__':
    # Parse options
    from optparse import OptionParser
    from phonopy.interface.vasp import read_vasp

    parser = OptionParser()
    parser.set_defaults( atom_num = 1,
                         pos_str = None )
    parser.add_option("-n", dest="atom_num", type="int",
			help="center atom")
    parser.add_option("-p", dest="pos_str", type="string",
			help="center position")
    (options, args) = parser.parse_args()


    cell = read_vasp('SPOSCAR')

    if options.pos_str == None:
        center = cell.get_scaled_positions()[ options.atom_num - 1 ]
    else:
        center = np.array( [ float(x) for x in options.pos_str.split() ] )

    for i in range(cell.get_number_of_atoms()):
        print "%2d  %f" % (i+1, get_distance_with_center(cell, center, i))
