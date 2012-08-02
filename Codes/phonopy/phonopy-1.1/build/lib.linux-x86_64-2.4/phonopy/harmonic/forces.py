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

class Forces:
    """
    forces: Forces on atoms in a supercell with a displacement in Cartesian coordinate
      [ [ F_1x, F_1y, F_1z ], 
        [ F_2x, F_2y, F_2z ], 
        ... ]
    displacement: An atomic displacement in Cartesian coordiante
      [ d_x, d_y, d_z ]
    """
    
    def __init__(self, atom_number, displacement, forces,
                 is_translational_invariance=False):
        self.atom_number = atom_number
        self.displacement = displacement
        self.forces = np.array(forces)
        if is_translational_invariance:
            self.set_translational_invariance()

    def get_atom_number(self):
        return self.atom_number

    def get_displacement(self):
        return self.displacement

    def get_forces(self):
        return self.forces

    def set_translational_invariance(self):
        self.forces = self.forces - \
            np.sum(self.forces, axis=0) / self.forces.shape[0]
        
