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

directions_axis = np.array([[ 1, 0, 0 ],
                            [ 0, 1, 0 ],
                            [ 0, 0, 1 ]])

directions_diag = np.array([[ 1, 0, 0 ],
                            [ 0, 1, 0 ],
                            [ 0, 0, 1 ],
                            [ 1, 1, 0 ],
                            [ 1, 0, 1 ],
                            [ 0, 1, 1 ],
                            [ 1,-1, 0 ],
                            [ 1, 0,-1 ],
                            [ 0, 1,-1 ],
                            [ 1, 1, 1 ],
                            [ 1, 1,-1 ],
                            [ 1,-1, 1 ],
                            [-1, 1, 1 ]])

def get_least_displacements( symmetry,
                             is_plusminus='auto',
                             is_diagonal=True,
                             log_level=0 ):
    """
    Return least displacements

    Format:
      [[ atom_num, disp_x, disp_y, disp_z ],
       [ atom_num, disp_x, disp_y, disp_z ],
       ...
      ]
    """
    displacements = []
    if is_diagonal:
        directions = directions_diag
    else:
        directions = directions_axis

    if log_level > 2:
        print "Site point symmetry:"

    for atom_num in symmetry.get_independent_atoms():
        site_symmetry = symmetry.get_site_symmetry(atom_num)

        if log_level > 2:
            print "Atom %d" % ( atom_num + 1 )
            for i, rot in enumerate( site_symmetry ):
                print "----%d----" % ( i + 1 )
                for v in rot:
                    print "%2d %2d %2d" % tuple(v)
        
        for disp in get_displacement(site_symmetry, directions):
            displacements.append([atom_num, disp[0], disp[1], disp[2]])
            if is_plusminus == 'auto':
                if is_minus_displacement(disp, site_symmetry):
                    displacements.append([atom_num, -disp[0], -disp[1], -disp[2]])
            elif is_plusminus == True:
                displacements.append([atom_num, -disp[0], -disp[1], -disp[2]])

    return displacements

def get_displacement(site_symmetry, directions=directions_diag):
    # Number of least displacements is,
    # One
    disp = get_displacement_one(site_symmetry, directions)
    if not disp == None:
        return disp
    # Two
    disp = get_displacement_two(site_symmetry, directions)
    if not disp == None:
        return disp
    # Three
    return [directions[0], directions[1], directions[2]]

def get_displacement_one(site_symmetry, directions=directions_diag):
    for direction in directions:
        rot_directions = []
        for rot in site_symmetry:
            rot_directions.append( np.dot( direction, rot.transpose() ) )
        num = site_symmetry.shape[0]
        for i in range(num):
            for j in range(i+1, num):
                for k in range(j+1, num):
                    det = np.linalg.det(np.array([rot_directions[i], rot_directions[j],
                                                  rot_directions[k]]))
                    if det != 0:
                        return [direction]
    return None

def get_displacement_two(site_symmetry, directions=directions_diag):
    for direction in directions:
        rot_directions = []
        for rot in site_symmetry:
            rot_directions.append( np.dot( direction, rot.transpose() ) )
        num = site_symmetry.shape[0]
        for i in range(num):
            for j in range(i+1,num):
                for sec_direction in directions:
                    det = np.linalg.det(np.array([rot_directions[i], rot_directions[j],
                                                  sec_direction]))
                    if det != 0:
                        return [direction, sec_direction]
    return None
    
def is_minus_displacement(direction, site_symmetry):
    is_minus = True
    for sym in site_symmetry:
        rot_direction = np.dot(direction, sym.transpose())
        if (rot_direction + direction).any():
            continue
        else:
            is_minus = False
            break
    return is_minus
        
def print_displacements(symmetry, directions=directions_diag):
    displacements = get_least_displacements(symmetry, directions)
    print "Least displacements:"
    print "Atom       Directions"
    print "----------------------------"
    for key in displacements:
        print "%4d  " % (key+1),
        for x in displacements[key]:
            print x,
        print

