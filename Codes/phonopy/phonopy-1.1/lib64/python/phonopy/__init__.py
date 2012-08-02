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
from phonopy.structure.atoms import Atoms
from phonopy.structure.symmetry import Symmetry
from phonopy.structure.cells import get_supercell, Primitive, print_cell
from phonopy.harmonic.displacement import get_least_displacements, print_displacements
from phonopy.harmonic.force_constant import get_force_constants, symmetrize_force_constants, rotational_invariance
from phonopy.harmonic.dynamical_matrix import DynamicalMatrix, DynamicalMatrixNAC
from phonopy.phonon.band_structure import BandStructure
from phonopy.phonon.thermal_properties import ThermalProperties
from phonopy.harmonic.forces import Forces
from phonopy.phonon.mesh import Mesh
from phonopy.units import VaspToTHz
from phonopy.phonon.dos import TotalDos, PartialDos
from phonopy.phonon.thermal_displacement import ThermalDisplacements, ThermalDistances
from phonopy.phonon.animation import Animation
from phonopy.phonon.modulation import write_modulations
from phonopy.phonon.qpoints_mode import write_yaml as write_yaml_qpoints

class Phonopy:
    def __init__( self,
                  unitcell,
                  supercell_matrix,
                  is_preprocess=True,
                  distance=0.01,
                  symprec=1e-5,
                  factor=VaspToTHz,
                  is_nosym=False,
                  log_level=0 ):
        self.symprec = symprec
        self.unitcell = unitcell
        self.supercell_matrix = supercell_matrix
        self.distance = distance
        self.factor = factor
        self.is_nosym = is_nosym
        self.log_level = log_level
        self.supercell = None
        self.__supercell()
        self.symmetry = None
        self.__symmetry()

        if is_preprocess:
            self.displacements = None
            self.displacement_directions = None
            self.supercells_with_displacements = None
            self.set_displacements()

        # set_post_process
        self.primitive = None
        self.dynamical_matrix = None
        self.is_nac = False

        # set_force_constants or set_forces
        self.force_constants = None

        # set_band_structure
        self.__band_structure = None

        # set_mesh
        self.__mesh = None

        # set_thermal_properties
        self.__thermal_properties = None

        # set_thermal_displacements
        self.__thermal_displacements = None

        # set_partial_DOS
        self.__pdos = None

        # set_total_DOS
        self.__total_dos = None
        
    def __supercell(self):
        self.supercell = get_supercell( self.unitcell,
                                        self.supercell_matrix,
                                        self.symprec )

    def get_supercell( self ):
        return self.supercell

    def set_supercell( self, supercell ):
        self.supercell = supercell

    def get_primitive( self ):
        return self.primitive

    def set_primitive( self, primitive ):
        self.primitive = primitive

    def __symmetry(self):
        self.symmetry = Symmetry( self.supercell,
                                  self.symprec,
                                  self.is_nosym )

    def get_symmetry(self):
        return self.symmetry

    def set_displacements( self,
                           is_plusminus='auto',
                           is_diagonal=True ):
        """
        displacements:
          List of displacements in Cartesian coordinates.
          See 'set_special_displacements'
        
        displacement_directions:
          List of directions with respect to axes. This gives only the
          symmetrically non equivalent directions. The format is like:
             [[ 0, 1, 0, 0 ],
              [ 7, 1, 0, 1 ], ...]
          where each list is defined by:
             First value:      Atom index in supercell starting with 0
             Second to fourth: If the direction is displaced or not ( 1, 0, or -1 )
                               with respect to the axes.
        """

        lattice = self.supercell.get_cell()
        self.displacements = []
        self.displacement_directions = \
            get_least_displacements( self.symmetry, 
                                     is_plusminus=is_plusminus,
                                     is_diagonal=is_diagonal,
                                     log_level=self.log_level )

        for disp in self.displacement_directions:
            atom_num = disp[0]
            disp_cartesian = np.dot(disp[1:], lattice)
            disp_cartesian *= self.distance / np.linalg.norm(disp_cartesian)
            self.displacements.append( [ atom_num,
                                         disp_cartesian[0],
                                         disp_cartesian[1],
                                         disp_cartesian[2] ] )

        self.__supercells_with_displacements()

    def set_special_displacements(self, displacements):
        """
        This method orverwrites displacements that were automatically
        determined in the post-process.

        displacemsts: List of disctionaries
           [[ 0, 0.01, 0.00, 0.00 ], ... ]
        where each set of elements is defined by:
           First value:      Atom index in supercell starting with 0
           Second to fourth: Displacement in Cartesian coordinates
        """

        self.displacements = displacements
        self.__supercells_with_displacements()

    def get_displacements( self ):
        return self.displacements

    def get_displacement_directions( self ):
        return self.displacement_directions

    def print_displacements(self):
        print "Least displacements:"
        print " Atom       Displacement"
        print " ----------------------------"
        for disp in self.displacements:
            print " %4d  " % disp[0],
            print disp[1:4]

    def __supercells_with_displacements(self):
        supercells = []
        for disp in self.displacements:
            positions = self.supercell.get_positions()
            positions[disp[0]] += disp[1:4]
            supercells.append( Atoms( 
                    numbers = self.supercell.get_atomic_numbers(),
                    masses = self.supercell.get_masses(),
                    positions = positions,
                    cell = self.supercell.get_cell(),
                    pbc = True ) )

        self.supercells_with_displacements = supercells
                        
    def get_supercells_with_displacements(self):
        return self.supercells_with_displacements

    def set_post_process( self,
                          primitive_matrix,
                          set_of_forces=None,
                          force_constants=None,
                          is_nac=False ):
        """
        Set forces to prepare phonon calculations. The order of
        'set_of_forces' has to correspond to that of 'displacements'.

        primitive_matrix:
          Relative axes of primitive cell to the input unit cell.
          Relative axes to the supercell is calculated by:
             supercell_matrix^-1 * primitive_matrix
          Therefore primitive cell lattice is finally calculated by:
             ( supercell_lattice * ( supercell_matrix )^-1 * primitive_matrix )^T

        set_of_forces:
           [ [ [ f_1x, f_1y, f_1z ], [ f_2x, f_2y, f_2z ], ... ], # first supercell
             [ [ f_1x, f_1y, f_1z ], [ f_2x, f_2y, f_2z ], ... ], # second supercell
             ...                                                   ]
        """

        self.is_nac = is_nac

        if not set_of_forces == None:
            self.set_forces( set_of_forces )
        elif not force_constants == None:
            self.set_force_constants( force_constants )
        elif self.force_constants == None:
            print "In set_post_process, set_of_forces or force_constants"
            print "has to be set."
            sys.exit(1)

        # Primitive cell
        inv_supercell_matrix = np.linalg.inv( self.supercell_matrix )
        self.primitive = Primitive( self.supercell,
                                    np.dot( inv_supercell_matrix, primitive_matrix ),
                                    self.symprec )

        # Dynamical Matrix
        if self.is_nac:
            self.dynamical_matrix = \
                DynamicalMatrixNAC( self.supercell,
                                    self.primitive,
                                    self.force_constants,
                                    symprec=self.symprec )
        else:
            self.dynamical_matrix = \
                DynamicalMatrix( self.supercell,
                                 self.primitive,
                                 self.force_constants,
                                 symprec=self.symprec )

    def set_nac_params( self, nac_params, method='wang' ):
        if self.is_nac:
            self.dynamical_matrix.set_nac_params( nac_params, method )

    def get_dynamical_matrix( self ):
        return self.dynamical_matrix

    def set_forces( self,
                    set_of_forces,
                    is_tensor_symmetry=False ):
        # Forces
        forces = []
        for i, disp in enumerate( self.displacements ):
            forces.append( Forces( disp[0],
                                   disp[1:4],
                                   set_of_forces[i] ) )

        # Force constants
        self.force_constants = get_force_constants( forces,
                                                    self.symmetry,
                                                    self.supercell,
                                                    is_tensor_symmetry )

    def set_force_constants( self, force_constants ):
        self.force_constants = force_constants

    def symmetrize_force_constants( self, iteration=3 ):
        symmetrize_force_constants( self.force_constants, iteration )
        
    def get_force_constants( self ):
        return self.force_constants

    def get_rotational_condition_of_fc( self ):
        return rotational_invariance( self.force_constants,
                                      self.supercell,
                                      self.primitive,
                                      self.symprec )

    # Frequency at a q-point
    def get_frequencies(self, q):
        """
        Calculate phonon frequency
        
        q: k-vector in reduced coordinates of primitive cell
        """
        self.dynamical_matrix.set_dynamical_matrix(q)
        dm = self.dynamical_matrix.get_dynamical_matrix()
        frequencies = []
        for eig in np.linalg.eigvalsh(dm):
            if eig < 0:
                frequencies.append(-np.sqrt(-eig))
            else:
                frequencies.append(np.sqrt(eig))
            
        return np.array(frequencies) * self.factor

    # Band structure
    def set_band_structure( self,
                            bands,
                            is_eigenvectors=False ):

        self.__band_structure = BandStructure( bands,
                                               self.dynamical_matrix,
                                               self.primitive,
                                               is_eigenvectors=is_eigenvectors,
                                               factor=self.factor )

    def get_band_structure( self ):
        band = self.__band_structure
        return ( band.get_distances(),
                 band.get_qpoints(),
                 band.get_eigenvalues(),
                 band.get_eigenvectors() )

    def plot_band_structure( self, symbols=None ):
        return self.__band_structure.plot_band( symbols )

    def write_yaml_band_structure( self ):
        self.__band_structure.write_yaml()

    # Mesh sampling
    def set_mesh( self,
                  mesh,
                  shift=None,
                  is_time_reversal=True,
                  is_symmetry=True,
                  is_eigenvectors=False ):

        self.__mesh = Mesh( self.dynamical_matrix,
                            self.primitive,
                            mesh,
                            shift=shift,
                            is_time_reversal=is_time_reversal,
                            is_symmetry=is_symmetry,
                            is_eigenvectors=is_eigenvectors,
                            factor=self.factor,
                            symprec=self.symprec )

    def get_mesh( self ):
        return ( self.__mesh.get_weights(),
                 self.__mesh.get_qpoints(),
                 self.__mesh.get_eigenvalues(),
                 self.__mesh.get_eigenvectors() )

    def write_yaml_mesh( self ):
        self.__mesh.write_yaml()

    # Thermal property
    def set_thermal_properties( self, t_step=10, t_max=1000, t_min=0,
                                cutoff_eigenvalue=None ):
        if self.__mesh==None:
            print "set_mesh has to be done before set_thermal_properties"
            sys.exit(1)

        tp = ThermalProperties( self.__mesh.get_eigenvalues(),
                                weights=self.__mesh.get_weights(),
                                factor=self.factor,
                                cutoff_eigenvalue=cutoff_eigenvalue )
        tp.set_thermal_properties( t_step, t_max, t_min )
        self.__thermal_properties = tp

    def get_thermal_properties( self ):
        temps, fe, entropy, cv = \
            self.__thermal_properties.get_thermal_properties()
        return np.array([ temps, fe, entropy, cv ]).transpose()

    def plot_thermal_properties( self ):
        return self.__thermal_properties.plot_thermal_properties()

    def write_yaml_thermal_properties( self ):
        self.__thermal_properties.write_yaml()

    # Partial DOS
    def set_partial_DOS( self,
                         sigma=None,
                         omega_min=None,
                         omega_max=None,
                         omega_pitch=None ):

        if self.__mesh==None:
            print "set_mesh has to be done before set_thermal_properties"
            sys.exit(1)
        if self.__mesh.get_eigenvectors() == None:
            print "Eigenvectors have to be calculated."
            sys.exit(1)

        pdos = PartialDos( self.__mesh.get_eigenvalues(),
                           self.__mesh.get_weights(),
                           self.__mesh.get_eigenvectors(),
                           factor=self.factor,
                           sigma=sigma )
        pdos.set_draw_area( omega_min,
                            omega_max,
                            omega_pitch )
        pdos.calculate()
        self.__pdos = pdos

    def get_partial_DOS( self ):
        """
        Retern omegas and partial_dos.
        The first element is omegas and the second is partial_dos.
        
        omegas: [ freq1, freq2, ... ]
        partial_dos:
          [[elem1-freq1, elem1-freq2, ... ],
           [elem2-freq1, elem2-freq2, ... ],
           ... ]

          where
           elem1: atom1-x compornent
           elem2: atom1-y compornent
           elem3: atom1-z compornent
           elem4: atom2-x compornent
           ...
        """
        return self.__pdos.get_partial_dos()

    def plot_partial_DOS( self, pdos_indices ):
        return self.__pdos.plot_pdos( pdos_indices )

    def write_partial_DOS( self ):
        self.__pdos.write()

    # Total DOS
    def set_total_DOS( self,
                       sigma=None,
                       omega_min=None,
                       omega_max=None,
                       omega_pitch=None ):

        if self.__mesh==None:
            print "set_mesh has to be done before set_thermal_properties"
            sys.exit(1)

        total_dos = TotalDos( self.__mesh.get_eigenvalues(),
                              self.__mesh.get_weights(),
                              factor=self.factor,
                              sigma=sigma )
        total_dos.set_draw_area( omega_min,
                                 omega_max,
                                 omega_pitch )
        total_dos.calculate()
        self.__total_dos = total_dos

    def get_total_DOS( self ):
        """
        Retern omegas and total dos.
        The first element is omegas and the second is total dos.
        
        omegas: [ freq1, freq2, ... ]
        total_dos: [ dos1, dos2, ... ]
        """
        return self.__total_dos.get_dos()

    def plot_total_DOS( self ):
        return self.__total_dos.plot_dos()

    def write_total_DOS( self ):
        self.__total_dos.write()

    # Thermal displacement
    def set_thermal_displacements( self,
                                   t_step=10,
                                   t_max=1000,
                                   t_min=0,
                                   projector=None,
                                   cutoff_eigenvalue=None ):
        """
        cutoff_eigenvalue:
          phonon modes that have frequencies below cutoff_eigenvalue
          are ignored.
          e.g. 0.1 (THz)
        """
        if self.__mesh==None:
            print "set_mesh has to be done before set_thermal_properties"
            sys.exit(1)
        if self.__mesh.get_eigenvectors() == None:
            print "Eigenvectors have to be calculated."
            sys.exit(1)

        td = ThermalDisplacements( self.__mesh.get_eigenvalues(),
                                   self.__mesh.get_eigenvectors(),
                                   self.__mesh.get_weights(),
                                   self.primitive.get_masses(),
                                   factor=self.factor,
                                   cutoff_eigenvalue=cutoff_eigenvalue )
        td.set_temperature_range( t_min, t_max, t_step )
        if not projector==None:
            td.project_eigenvectors( projector, self.primitive.get_cell() )
        td.set_thermal_displacements()
        
        self.__thermal_displacements = td

    def plot_thermal_displacements( self, is_legend=False ):
        return self.__thermal_displacements.plot_thermal_displacements( is_legend )

    def write_yaml_thermal_displacements( self ):
        self.__thermal_displacements.write_yaml()

    # Thermal displacement
    def set_thermal_distances( self,
                               atom_pairs,
                               t_step=10,
                               t_max=1000,
                               t_min=0,
                               cutoff_eigenvalue=None ):
        """
        atom_pairs: List of list
          Mean square distances are calculated for the atom_pairs
          e.g. [ [ 1, 2 ], [ 1, 4 ] ]

        cutoff_eigenvalue:
          phonon modes that have frequencies below cutoff_eigenvalue
          are ignored.
          e.g. 0.1 (THz)
        """

        td = ThermalDistances( self.__mesh.get_eigenvalues(),
                               self.__mesh.get_eigenvectors(),
                               self.__mesh.get_weights(),
                               self.supercell,
                               self.primitive,
                               self.__mesh.get_qpoints(),
                               symprec=self.symprec,
                               factor=self.factor,
                               cutoff_eigenvalue=cutoff_eigenvalue )
        td.set_temperature_range( t_min, t_max, t_step )
        td.set_thermal_distances( atom_pairs )

        self.__thermal_distances = td

    def write_yaml_thermal_distances( self ):
        self.__thermal_distances.write_yaml()


    # Q-points mode
    def write_yaml_qpoints( self,
                            qpoints,
                            is_eigenvectors=False,
                            factor=VaspToTHz ):
        
        write_yaml_qpoints( qpoints,
                            self.primitive,
                            self.dynamical_matrix,
                            is_eigenvectors=is_eigenvectors,
                            factor=self.factor )

    # Animation
    def write_animation( self,
                         qpoint=None,
                         anime_type='v_sim',
                         band_index=None,
                         amplitude=None,
                         num_div=None,
                         shift=None ):
        if qpoint==None:
            animation = Animation( [0, 0, 0],
                                   self.dynamical_matrix,
                                   self.primitive,
                                   shift=shift )
        else:
            animation = Animation( qpoint,
                                   self.dynamical_matrix,
                                   self.primitive,
                                   shift=shift )
        if anime_type=='v_sim':
            animation.write_v_sim( self.factor )

            
        if ( anime_type=='arc' or
             anime_type=='xyz' or
             anime_type=='jmol' or
             anime_type=='poscar' ):
            if band_index==None or amplitude==None or num_div==None:
                print "Parameters are not correctly set for animation."
                sys.exit(1)

            if anime_type=='arc' or anime_type==None:
                animation.write_arc( band_index,
                                     amplitude,
                                     num_div )
    
            if anime_type=='xyz':
                animation.write_xyz( band_index,
                                     amplitude,
                                     num_div,
                                     self.factor )
    
            if anime_type=='jmol':
                animation.write_xyz_jmol( amplitude=amplitude,
                                          factor=self.factor )
    
            if anime_type=='poscar':
                animation.write_POSCAR( band_index,
                                        amplitude,
                                        num_div )


    def write_modulation( self, setting ):
        write_modulations( self.dynamical_matrix,
                           self.primitive,
                           setting )
