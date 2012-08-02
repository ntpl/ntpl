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

import sys
import numpy as np
import phonopy.structure.spglib as spg
from phonopy.structure.atoms import Atoms

def find_primitive( cell, symprec=1e-5 ):
    """
    A primitive cell is searched in the input cell. When a primitive
    cell is found, an object of Atoms class of the primitive cell is
    returned. When not, None is returned.
    """
    lattice, positions, numbers = spg.find_primitive( cell, symprec )
    if lattice == None:
        return None
    else:
        return Atoms( numbers=numbers,
                      scaled_positions=positions,
                      cell=lattice,
                      pbc=True )

class Symmetry:
    def __init__(self, cell, symprec=1e-5, is_nosym=False):
        self.__cell = cell
        self.symprec = symprec

        self.symmetry_operations = None
        self.international_table = None
        self.dataset = None
        self.wyckoff_letters = None
        self.__map_atoms = None
        if is_nosym:
            self.__set_nosym()
        else:
            self.__symmetry_dataset()

        self.pointgroup_operations = None
        self.__pointgroup_operations()

        self.independent_atoms = None
        self.map_operations = None
        self.__mapping_atoms()

    def get_symmetry_operations(self):
        return self.symmetry_operations

    def get_symmetry_operation(self, operation_number):
        operation = self.symmetry_operations
        return {'rotations': operation['rotations'][operation_number],
                'translations': operation['translations'][operation_number]}

    def __pointgroup_operations(self):
        rotations = []
        for rot in self.symmetry_operations['rotations']:
            is_same = False
            for tmp_rot in rotations:
                if (tmp_rot==rot).all():
                    is_same = True
                    break
            if not is_same:
                rotations.append(rot)

        self.pointgroup_operations = np.array(rotations)

    def get_pointgroup_operations(self):
        return self.pointgroup_operations

    def __symmetry_dataset(self):
        self.dataset = spg.get_symmetry_dataset( self.__cell, self.symprec )
        self.symmetry_operations = \
            { 'rotations': self.dataset['rotations'],
              'translations': self.dataset['translations'] }
        self.international_table = "%s (%d)" % ( self.dataset['international'],
                                                 self.dataset['number'] )
        self.wyckoff_letters = self.dataset['wyckoffs']
        self.map_atoms = self.dataset['equivalent_atoms']

    def get_international_table(self):
        return self.international_table

    def get_Wyckoff_letters(self):
        return self.wyckoff_letters

    def get_dataset( self ):
        """
        Detail of dataset is found in spglib.get_symmetry_dataset.
        """
        return self.dataset

    def __mapping_atoms(self):
        ops = self.symmetry_operations
        pos = self.__cell.get_scaled_positions()
        map_operations = np.zeros( len( pos ), dtype=int )
        independent_atoms = []

        for i, eq_atom in enumerate( self.map_atoms ):
            if i==eq_atom:
                independent_atoms.append(i)
            for j, (r, t) in enumerate(
                zip( ops['rotations'], ops['translations'] ) ):
                
                diff = np.dot( pos[i], r.T ) + t - pos[eq_atom]
                if ( abs( diff - diff.round() ) < self.symprec ).all():
                    map_operations[ i ] = j
                    break

        self.independent_atoms = np.array( independent_atoms )
        self.map_operations = map_operations

    def get_independent_atoms(self):
        return self.independent_atoms

    def get_map_atoms(self):
        return self.map_atoms

    def get_map_operations(self):
        return self.map_operations

    def get_site_symmetry(self, atom_number):
        pos = self.__cell.get_scaled_positions()[atom_number]
        symprec = self.symprec
        rot = self.symmetry_operations['rotations']
        trans = self.symmetry_operations['translations']
        site_symmetries = []

        for r, t in zip( rot, trans ):
            rot_pos = np.dot( pos, r.T ) + t
            diff = pos - rot_pos
            if ( abs( diff - diff.round() ) < symprec ).all():
                site_symmetries.append( r )

        return np.array(site_symmetries)

    def get_symmetry_tolerance(self):
        return self.symprec

    def __set_nosym(self):
        translations = []
        rotations = []
        s2u_map = self.__cell.get_supercell_to_unitcell_map()
        positions = self.__cell.get_scaled_positions()

        for i, j in enumerate( s2u_map ):
            if j==0:
                ipos0 = i
                break

        for i, p in zip( s2u_map, positions ):
            if i==0:
                trans = p - positions[ipos0]
                trans -= np.floor( trans )
                translations.append( trans )
                rotations.append( np.eye( 3, dtype=int ) )

        self.symmetry_operations = { 'rotations': np.array( rotations ),
                                     'translations': np.array( translations ) }
        self.international_table = 'P1 (1)'
        self.wyckoff_letters = ['a'] * self.__cell.get_number_of_atoms()
        self.map_atoms = s2u_map


def get_ir_reciprocal_mesh(mesh,
                           cell,
                           is_shift=np.zeros(3, dtype=int),
                           is_time_reversal=False,
                           symprec=1e-5):
    """
    Return k-point map to the irreducible k-points and k-point grid points .
    The symmetry is serched from the input cell.
    is_shift=[ 0, 0, 0 ] gives Gamma center mesh.
    """

    return spg.get_ir_reciprocal_mesh(mesh,
                                      cell,
                                      is_shift,
                                      is_time_reversal,
                                      symprec)

def get_ir_kpoints(kpoints, cell, is_time_reversal=False, symprec=1e-5):

    return spg.get_ir_kpoints(kpoints, cell, is_time_reversal, symprec)

    


if __name__ == '__main__':
    # Parse options
    import sys
    from optparse import OptionParser
    from phonopy.interface.vasp import read_vasp

    class TestSymmetry( Symmetry ):
        def __init__( self, cell ):
            Symmetry.__init__( self, cell )
            
        def check_consistency( self ):
            count = 0
            is_consistent = True
            for r1, t1 in zip( self.dataset['rotations'],
                               self.dataset['translations'] ):
                for r2, t2 in zip( self.symmetry_operations['rotations'],
                                   self.symmetry_operations['translations'] ):
                    if ( (r1 - r2) == 0 ).all():
                        diff = t1 - t2
                        diff -= diff.round()
                        if ( abs(diff) < 1e-5 ).all():
                            count += 1
                            break
            if not ( count== len( self.symmetry_operations['rotations'] ) ):
                print "Symmetry operations may not organize a group."
                is_consistent = False
    
            return is_consistent


    cell = read_vasp(sys.argv[1])
    symmetry = TestSymmetry( cell )

    print symmetry.check_consistency()
    print symmetry.get_dataset()['transformation_matrix']
    print symmetry.get_dataset()['origin_shift']
