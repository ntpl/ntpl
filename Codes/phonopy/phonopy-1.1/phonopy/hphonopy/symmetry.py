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
from phonopy.structure.symmetry import Symmetry, find_primitive
from phonopy.structure.cells import Primitive, print_cell, get_supercell
from phonopy.interface.vasp import write_vasp
from phonopy.structure.atoms import Atoms

def get_symmetry_yaml( cell, symmetry, phonopy_version=None ):
    rotations = symmetry.get_symmetry_operations()['rotations']
    translations = symmetry.get_symmetry_operations()['translations']
    
    atom_sets = symmetry.get_map_atoms()
    independent_atoms = symmetry.get_independent_atoms()
    wyckoffs = symmetry.get_Wyckoff_letters()

    yaml = ""

    if not phonopy_version==None:
        yaml += "phonopy_version: %s\n" % phonopy_version

    yaml += "space_group_type: " + symmetry.get_international_table() + "\n"
    yaml += "space_group_operations:\n"
    for i, (r, t) in enumerate( zip( rotations, translations ) ):
        yaml += "- rotation: # %d\n" % (i+1)
        for vec in r:
            yaml += "  - [%2d%3d%3d]\n" % tuple( vec )
        yaml += "  translation: [%8.5f,%8.5f,%8.5f]\n" % tuple(t)

    yaml += "atom_mapping:\n"
    for i, atom_num in enumerate( atom_sets ):
        yaml += "  %d: %d\n" % ( i+1, atom_num+1 )
    yaml += "site_symmetries:\n"
    for i in independent_atoms:
        yaml += "- atom: %d\n" % (i+1)
        yaml += "  Wyckoff: %s\n" % ( wyckoffs[i] )
        yaml += "  rotations:\n"
        for j, r in enumerate( symmetry.get_site_symmetry( i ) ):
            yaml += "  - # %d\n" % (j+1)
            for vec in r:
                yaml += "    - [%2d%3d%3d]\n" % tuple( vec )

    return yaml

def check_symmetry( input_cell,
                    primitive_axis=np.eye( 3, dtype=float ),
                    symprec=1e-5,
                    phonopy_version=None ):

    cell = Primitive( input_cell,
                      primitive_axis,
                      symprec )

    symmetry = Symmetry( cell, symprec )
    print get_symmetry_yaml( cell, symmetry, phonopy_version ),

    primitive = find_primitive( cell, symprec )
    if not primitive==None:
        print "# Primitive cell was found. It is written into PPOSCAR."
        write_vasp( 'PPOSCAR', primitive )
        
        # Overwrite symmetry and cell
        symmetry = Symmetry( primitive, symprec )
        cell = primitive

    bravais_lattice, bravais_pos, bravais_numbers = \
        spg.refine_cell( cell, symprec )
    bravais = Atoms( numbers=bravais_numbers,
                     scaled_positions=bravais_pos,
                     cell=bravais_lattice,
                     pbc=True )
    print "# Bravais lattice is written into BPOSCAR."
    write_vasp( 'BPOSCAR', bravais )

