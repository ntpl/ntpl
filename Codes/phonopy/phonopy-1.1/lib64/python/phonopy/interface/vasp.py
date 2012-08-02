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
from phonopy.structure.atoms import Atoms, symbol_map, atom_data

class VasprunWrapper(object):
    """VasprunWrapper class
    This is used to avoid VASP 5.2.8 vasprun.xml defect at PRECFOCK,
    xml parser stops with error.
    """
    def __init__(self, filename):
        self.f = open(filename)

    def read(self, size=None):
        element = self.f.next()
        if element.find("PRECFOCK") == -1:
            return element
        else:
            return "<i type=\"string\" name=\"PRECFOCK\"></i>"

def write_supercells_with_displacements( supercell,
                                         cells_with_displacements ):
    write_vasp("SPOSCAR", supercell, direct=True)
    for i, cell in enumerate( cells_with_displacements ):
        write_vasp('POSCAR-%03d' % (i+1), cell, direct=True)

def get_forces_vasprun_xml( vasprun ):
    """
    vasprun = etree.iterparse( filename, tag='varray' )
    """
    forces = []
    num_atom = 0
    for event, element in vasprun:
        if element.attrib['name'] == 'forces':
            for v in element.xpath('./v'):
                forces.append([ float(x) for x in v.text.split() ])
                
    return np.array( forces )

def get_force_constants_vasprun_xml( vasprun ):
    fc_tmp = None
    num_atom = 0
    for event, element in vasprun:
        if num_atom==0:
            atom_types, masses, num_atom = get_atom_types_from_vasprun_xml( element )

        # Get Hessian matrix (normalized by masses)
        if element.tag == 'varray':
            if element.attrib['name'] == 'hessian':
                fc_tmp = []
                for v in element.xpath('./v'):
                    fc_tmp.append([ float(x) for x in v.text.strip().split() ])

    if fc_tmp==None:
        return False
    else:
        fc_tmp = np.array(fc_tmp)
        if not fc_tmp.shape==( num_atom*3, num_atom*3 ):
            return False
        # num_atom = fc_tmp.shape[0] / 3
        force_constants = np.zeros((num_atom, num_atom, 3, 3), dtype=float)
    
        for i in range(num_atom):
            for j in range(num_atom):
                force_constants[i, j] = fc_tmp[i*3:(i+1)*3, j*3:(j+1)*3]
    
        # Inverse normalization by atomic weights
        for i in range( num_atom ):
            for j in range( num_atom ):
                force_constants[i, j] *= -np.sqrt( masses[i] * masses[j] )
    
        return force_constants, atom_types

def get_atom_types_from_vasprun_xml( element ):
    atom_types = []
    masses = []
    num_atom = 0
    
    if element.tag == 'array':
        if 'name' in element.attrib:
            if element.attrib['name'] == 'atomtypes':
                for rc in element.xpath('./set/rc'):
                    atom_info = [ x.text for x in rc.xpath('./c') ]
                    num_atom += int( atom_info[0] )
                    atom_types.append( atom_info[1].strip() )
                    masses += ( [float( atom_info[2] )] * int( atom_info[0] ) )

    return atom_types, masses, num_atom

def get_force_constants_OUTCAR( filename ):
    file = open( filename )
    while 1:
        line = file.readline()
        if line == '':
            print "Force constants could not be found."
            return 0

        if line[:19] == " SECOND DERIVATIVES":
            break

    file.readline()
    num_atom = int( ( ( file.readline().split() )[-1].strip() )[:-1] )

    fc_tmp = []
    for i in range( num_atom * 3 ):
        fc_tmp.append( [ float(x) for x in ( file.readline().split() )[1:] ] )

    fc_tmp = np.array( fc_tmp )

    force_constants = np.zeros( (num_atom, num_atom, 3, 3), dtype=float )
    for i in range( num_atom ):
        for j in range( num_atom ):
            force_constants[i, j] = -fc_tmp[i*3:(i+1)*3, j*3:(j+1)*3]

    return force_constants
    



#
# read VASP POSCAR
#
def expand_symbols( num_atoms, symbols=None ):
    expanded_symbols = []
    is_symbols = True
    if symbols == None:
        is_symbols = False
    else:
        if not len( symbols ) == len( num_atoms ):
            is_symbols = False
        else:
            for s in symbols:
                if not s in symbol_map:
                    is_symbols = False
                    break
    
    if is_symbols:
        for s, num in zip( symbols, num_atoms ):
            expanded_symbols += [s] * num
    else:
        for i, num in enumerate( num_atoms ):
            expanded_symbols += [ atom_data[i+1][1] ] * num

    return expanded_symbols

def is_exist_symbols( symbols ):
    for s in symbols:
        if not ( s in symbol_map ):
            return False
    return True

def read_vasp( filename, symbols=None ):
    file = open( filename )
    
    lines = file.readlines()

    line1 = [ x for x in lines[0].split() ]
    if is_exist_symbols( line1 ):
        symbols = line1

    scale = float(lines[1])

    cell = []
    for i in range( 2, 5 ):
        cell.append( [ float(x) for x in lines[i].split()[:3] ] )
    cell = np.array( cell ) * scale

    try:
        num_atoms = np.array([ int(x) for x in lines[5].split() ])
        line_at = 6
    except ValueError:
        symbols = [ x for x in lines[5].split() ]
        num_atoms = np.array([ int(x) for x in lines[6].split() ])
        line_at = 7
    
    expaned_symbols = expand_symbols( num_atoms, symbols )

    if lines[ line_at ][0].lower() == 's':
        line_at += 1

    is_scaled = True
    if ( lines[ line_at ][0].lower() == 'c' or
         lines[ line_at ][0].lower() == 'k' ):
        is_scaled = False

    line_at += 1

    positions = []
    for i in range( line_at, line_at + num_atoms.sum() ):
        positions.append( [ float(x) for x in lines[i].split()[:3] ] )

    if is_scaled:
        atoms = Atoms( symbols=expaned_symbols,
                       cell=cell,
                       scaled_positions=positions )
    else:
        atoms = Atoms( symbols=expaned_symbols,
                       cell=cell,
                       positions=positions )
        
    return atoms
                   
#
# write vasp POSCAR
#
def get_reduced_symbols( symbols ):
    reduced_symbols = []
    for s in symbols:
        if not ( s in reduced_symbols ):
            reduced_symbols.append(s)
    return reduced_symbols

def sort_positions_by_symbols( symbols, positions ):
    reduced_symbols = get_reduced_symbols( symbols )
    sorted_positions = []
    num_atoms = np.zeros( len( reduced_symbols ), dtype=int )
    for i, rs in enumerate( reduced_symbols ):
        for s, p in zip( symbols, positions ):
            if rs==s:
                sorted_positions.append( p )
                num_atoms[i] += 1
    return num_atoms, reduced_symbols, np.array( sorted_positions )

def write_vasp( filename, atoms, direct=True ):

    num_atoms, symbols, scaled_positions = \
        sort_positions_by_symbols( atoms.get_chemical_symbols(),
                                   atoms.get_scaled_positions() )

    lines = ""     
    for s in symbols:
        lines += "%s " % s
    lines += "\n"
    lines += "   1.0\n"
    for a in atoms.get_cell():
        lines += " %22.16f%22.16f%22.16f\n" % tuple( a )
    lines += ("%4d" * len( num_atoms )) % tuple( num_atoms )
    lines += "\n"
    lines += "Direct\n"
    for vec in scaled_positions:
        for x in ( vec - vec.round() ):
            if float('%20.16f' % x) < 0.0:
                lines += "%20.16f" % ( x + 1.0 )
            else:
                lines += "%20.16f" % (x)
        lines += "\n"

    f = open( filename, 'w' )
    f.write( lines )

if __name__ == '__main__':
    import sys
    atoms = read_vasp( sys.argv[1] )
    write_vasp( '%s-new' % sys.argv[1], atoms )
