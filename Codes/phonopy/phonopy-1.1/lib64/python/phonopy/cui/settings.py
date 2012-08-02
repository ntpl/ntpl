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

def fracval( frac ):
    if frac.find('/') == -1:
        return float( frac )
    else:
        x = frac.split('/')
        return float( x[0] ) / float( x[1] )

class Settings:
    def __init__( self ):

        self.chemical_symbols = None
        self.is_eigenvectors = False
        self.is_diagonal_displacement = True
        self.is_plusminus_displacement = 'auto'
        self.is_tensor_symmetry = False
        self.is_translational_invariance = False
        self.is_rotational_invariance = False
        self.fc_symmetry_iteration = 0
        self.masses = None
        self.mesh = None
        self.omega_step = None
        self.primitive_matrix = np.eye(3, dtype=float)
        self.run_mode = None
        self.sigma = None
        self.supercell_matrix = None
        self.tmax = 1000
        self.tmin = 0
        self.tstep = 10

    def set_run_mode(self, run_mode):
        self.run_mode = run_mode

    def get_run_mode(self):
        return self.run_mode

    def set_supercell_matrix(self, matrix):
        self.supercell_matrix = matrix

    def get_supercell_matrix(self):
        return self.supercell_matrix

    def set_is_diagonal_displacement(self, is_diag):
        self.is_diagonal_displacement = is_diag

    def get_is_diagonal_displacement(self):
        return self.is_diagonal_displacement

    def set_is_plusminus_displacement(self, is_pm):
        self.is_plusminus_displacement = is_pm

    def get_is_plusminus_displacement(self):
        return self.is_plusminus_displacement

    def set_masses(self, masses):
        self.masses = masses

    def get_masses(self):
        return self.masses

    def set_chemical_symbols(self, symbols):
        self.chemical_symbols = symbols

    def get_chemical_symbols(self):
        return self.chemical_symbols

    def set_mesh_numbers(self, mesh):
        self.mesh = mesh

    def get_mesh_numbers(self):
        return self.mesh

    def set_mesh_shift(self, mesh_shift):
        self.mesh_shift = mesh_shift

    def get_mesh_shift(self):
        return self.mesh_shift

    def set_mesh_symmetry(self, mesh_symmetry=True):
        self.mesh_symmetry = mesh_symmetry

    def get_mesh_symmetry(self):
        return self.mesh_symmetry

    def set_time_symmetry(self, time_symmetry=True):
        self.time_symmetry = time_symmetry

    def get_time_symmetry(self):
        return self.time_symmetry

    def set_primitive_matrix(self, primitive_matrix):
        self.primitive_matrix = primitive_matrix

    def get_primitive_matrix(self):
        return self.primitive_matrix
        
    def set_is_eigenvectors(self, is_eigenvectors):
        self.is_eigenvectors = is_eigenvectors

    def get_is_eigenvectors(self):
        return self.is_eigenvectors

    def set_is_tensor_symmetry(self, is_tensor_symmetry):
        self.is_tensor_symmetry = is_tensor_symmetry

    def get_is_tensor_symmetry(self):
        return self.is_tensor_symmetry

    def set_is_translational_invariance(self, is_translational_invariance):
        self.is_translational_invariance = is_translational_invariance

    def get_is_translational_invariance(self):
        return self.is_translational_invariance

    def set_is_rotational_invariance(self, is_rotational_invariance):
        self.is_rotational_invariance = is_rotational_invariance

    def get_is_rotational_invariance(self):
        return self.is_rotational_invariance

    def set_fc_symmetry_iteration(self, iteration):
        self.fc_symmetry_iteration = iteration

    def get_fc_symmetry_iteration(self):
        return self.fc_symmetry_iteration

    def set_sigma(self, sigma):
        self.sigma = sigma

    def get_sigma(self):
        return self.sigma

    def set_omega_step(self, omega_step):
        self.omega_step = omega_step

    def get_omega_step(self):
        return self.omega_step

    def set_max_temperature( self, tmax ):
        self.tmax = tmax

    def get_max_temperature( self ):
        return self.tmax

    def set_min_temperature( self, tmin ):
        self.tmin = tmin

    def get_min_temperature( self ):
        return self.tmin

    def set_temperature_step( self, tstep ):
        self.tstep = tstep

    def get_temperature_step( self ):
        return self.tstep



# Parse phonopy setting filen
class ConfParser:
    def __init__(self, filename=None, options=None, option_list=None ):
        self.confs = {}
        self.tags = {}
        self.options = options
        self.option_list = option_list

        if not filename==None:
            self.__read_file( filename )
        if ( not options==None ) and ( not option_list==None ):
            self.__read_options()
        self.__parse_conf()

    def get_settings( self ):
        return self.settings

    def setting_error( self, message ):
        print message
        print "Please check the setting tags and options."
        sys.exit(1)

    def set_settings( self ):
        tags = self.tags

        # Is getting least displacements?
        if tags.has_key('create_displacements'):
            if tags['create_displacements'] == '.true.':
                self.settings.set_run_mode('displacements')
    
        # Supercell size
        if tags.has_key('supercell_matrix'):
            self.settings.set_supercell_matrix(tags['supercell_matrix'])

        # Atomic mass
        if tags.has_key('mass'):
            self.settings.set_masses(tags['mass'])
    
        # Chemical symbols
        if tags.has_key('atom_name'):
            self.settings.set_chemical_symbols(tags['atom_name'])
            
        # Diagonal displacement
        if tags.has_key('diag'):
            if tags['diag'] == '.false.':
                self.settings.set_is_diagonal_displacement(False)
    
        # Plus minus displacement
        if tags.has_key('pm'):
            if tags['pm'] == '.false.':
                self.settings.set_is_plusminus_displacement(False)
            if tags['pm'] == '.true.':
                self.settings.set_is_plusminus_displacement(True)
    
        # Primitive cell shape
        if tags.has_key('primitive_axis'):
            self.settings.set_primitive_matrix( tags['primitive_axis'] )
    
        # Is getting eigenvectors, too ?
        if tags.has_key('eigenvectors'):
            if tags['eigenvectors'] == '.true.':
                self.settings.set_is_eigenvectors(True)
    
        # Enforce force constant symmetry?
        if tags.has_key('fc_symmetry'):
            self.settings.set_fc_symmetry_iteration( int( tags['fc_symmetry'] ) )
    
        # Is translational invariance ?
        if tags.has_key('translation'):
            if tags['translation'] == '.true.':
                self.settings.set_is_translational_invariance(True)
    
        # Is rotational invariance ?
        if tags.has_key('rotational'):
            if tags['rotational'] == '.true.':
                self.settings.set_is_rotational_invariance(True)
    
        # Is force constants symmetry forced?
        if tags.has_key('tensor_symmetry'):
            if tags['tensor_symmetry'] == '.true.':
                self.settings.set_is_tensor_symmetry(True)

        # Mesh sampling numbers
        if tags.has_key('mesh_numbers'):
            self.settings.set_mesh_numbers( tags['mesh_numbers'] )

        # Spectram drawing step
        if tags.has_key('omega_step'):
            self.settings.set_omega_step(tags['omega_step'])

        # Smearing width
        if tags.has_key('sigma'):
            self.settings.set_sigma(tags['sigma'])

        # Temerature range
        if tags.has_key('tmax'):
            self.settings.set_max_temperature( tags['tmax'] )
        if tags.has_key('tmin'):
            self.settings.set_min_temperature( tags['tmin'] )
        if tags.has_key('tstep'):
            self.settings.set_temperature_step( tags['tstep'] )

    def __read_file(self, filename ):
        file = open( filename, 'r' )
        confs = self.confs
        is_continue = False
        for line in file:
            if line.strip() == '':
                is_continue = False
                continue
            
            if line.strip()[0] == '#':
                is_continue = False
                continue

            if is_continue:
                confs[left] += line.strip()
                confs[left] = confs[left].replace('+++', ' ')
                is_continue = False
                
            if line.find('=') != -1:
                left, right = [ x.strip().lower() for x in line.split('=') ]
                confs[left] = right

            if line.find('+++') != -1:
                is_continue = True

    def __read_options(self):
        for opt in self.option_list:

            if opt.dest=='is_displacement':
                if self.options.is_displacement:
                    self.confs['create_displacements'] = '.true.'

            if opt.dest=='is_plusminus_displacements':
                if self.options.is_plusminus_displacements:
                    self.confs['pm'] = '.true.'

            if opt.dest=='mesh_numbers':
                if not self.options.mesh_numbers==None:
                    self.confs['mesh_numbers'] = self.options.mesh_numbers

            if opt.dest=='primitive_axis':
                if not self.options.primitive_axis==None:
                    self.confs['primitive_axis'] = self.options.primitive_axis
                    
            if opt.dest=='supercell_dimension':
                if not self.options.supercell_dimension==None:
                    self.confs['dim'] = self.options.supercell_dimension

            if opt.dest=='is_nodiag':
                if self.options.is_nodiag:
                    self.confs['diag'] = '.false.'

            if opt.dest=='omega_step':
                if not self.options.omega_step==None:
                    self.confs['omega_step'] = self.options.omega_step

            if opt.dest=='sigma':
                if not self.options.sigma==None:
                    self.confs['sigma'] = self.options.sigma

            if opt.dest=='tmin':
                if not self.options.tmin==None:
                    self.confs['tmin'] = self.options.tmin

            if opt.dest=='tmax':
                if not self.options.tmax==None:
                    self.confs['tmax'] = self.options.tmax
                    
            if opt.dest=='tstep':
                if not self.options.tstep==None:
                    self.confs['tstep'] = self.options.tstep

    def __parse_conf(self):
        confs = self.confs

        for tag in confs.keys():
            if tag == 'create_displacements':
                self.set_parameter( 'create_displacements',
                                    confs['create_displacements'] )

            if tag == 'dim':
                matrix = [ int(x) for x in confs['dim'].split() ]
                if len( matrix ) == 9:
                    matrix = np.array( matrix ).reshape( 3, 3 )
                elif len( matrix ) == 3:
                    matrix = np.diag( matrix )
                else:
                    self.setting_error("Number of elements of DIM tag has to be 3 or 9.")

                if matrix.shape == ( 3, 3 ):
                    if np.linalg.det( matrix ) < 1:
                        self.setting_error('Determinant of supercell matrix has to be positive.')
                    else:
                        self.set_parameter( 'supercell_matrix', matrix )

            if tag == 'primitive_axis':
                if not len(confs['primitive_axis'].split()) == 9:
                    self.setting_error("Number of elements in PRIMITIVE_AXIS has to be 9.")
                p_axis = []
                for x in confs['primitive_axis'].split():
                    p_axis.append( fracval( x ) )
                p_axis = np.array( p_axis ).reshape(3,3)
                if np.linalg.det( p_axis ) < 1e-8:
                    self.setting_error("PRIMITIVE_AXIS has to have positive determinant.")
                self.set_parameter('primitive_axis', p_axis )

            if tag == 'mass':
                self.set_parameter('mass',
                                   [ float(x) for x in confs['mass'].split()])

            if tag == 'atom_name':
                self.set_parameter('atom_name',
                                   [ x.capitalize() for x in confs['atom_name'].split() ])
                        
            if tag == 'eigenvectors':
                self.set_parameter('eigenvectors',
                                     confs['eigenvectors'])

            if tag == 'diag':
                self.set_parameter('diag', confs['diag'])

            if tag == 'pm':
                self.set_parameter('pm', confs['pm'])

            if tag == 'translation':
                self.set_parameter('translation',
                                     confs['translation'])

            if tag == 'rotational':
                self.set_parameter('rotational',
                                     confs['rotational'])

            if tag == 'fc_symmetry':
                self.set_parameter('fc_symmetry',
                                     confs['fc_symmetry'])

            if tag == 'tensor_symmetry':
                self.set_parameter('tensor_symmetry',
                                     confs['tensor_symmetry'])

            if tag == 'mesh_numbers':
                vals = [ int(x) for x in confs['mesh_numbers'].split() ]
                if len(vals) < 3:
                    self.setting_error("Mesh numbers are incorrectly set.")
                self.set_parameter('mesh_numbers', vals[:3])

            if tag == 'omega_step':
                if isinstance( confs['omega_step'], str ):
                    val = float(confs['omega_step'].split()[0])
                else:
                    val = confs['omega_step']
                self.set_parameter('omega_step', val)

            if tag == 'sigma':
                if isinstance( confs['sigma'], str ):
                    val = float(confs['sigma'].split()[0])
                else:
                    val = confs['sigma']
                self.set_parameter('sigma', val)

            if tag == 'tmin':
                val = float(confs['tmin'].split()[0])
                self.set_parameter('tmin', val)

            if tag == 'tmax':
                val = float(confs['tmax'].split()[0])
                self.set_parameter('tmax', val)

            if tag == 'tstep':
                val = float(confs['tstep'].split()[0])
                self.set_parameter('tstep', val)

    def set_parameter( self, key, val ):
        self.tags[key] = val



#
# For phonopy
#
class PhonopySettings( Settings ):
    def __init__( self ):
        Settings.__init__( self )

        self.anime_band_index = None
        self.anime_amplitude = None
        self.anime_division = None
        self.anime_qpoint = None
        self.anime_shift = None
        self.anime_type = 'v_sim'
        self.dos = None
        self.dos_range = { 'min':  None,
                           'max':  None }
        self.thermal_atom_pairs = None
        self.is_dos_mode = False
        self.is_force_constants = False
        self.is_plusminus_displacement = 'auto'
        self.is_thermal_displacements = False
        self.is_thermal_distances = False
        self.is_thermal_properties = False
        self.modulation = None
        self.pdos_indices = None

    def set_run_mode(self, run_mode):
        modes = [ 'qpoints',
                  'mesh',
                  'band',
                  'anime',
                  'modulation',
                  'displacements' ]
        for mode in modes:
            if run_mode.lower() == mode:
                self.run_mode = run_mode

    def get_run_mode(self):
        return self.run_mode

    def set_is_force_constants(self, is_force_constants):
        self.is_force_constants = is_force_constants

    def get_is_force_constants(self):
        return self.is_force_constants

    def set_bands(self, bands):
        self.bands = bands

    def get_bands(self):
        return self.bands

    def set_mesh( self,
                  mesh,
                  mesh_shift=[0.,0.,0.],
                  time_symmetry=True,
                  mesh_symmetry=True ):
        self.mesh = mesh
        self.mesh_shift = mesh_shift
        self.time_symmetry = time_symmetry
        self.mesh_symmetry = mesh_symmetry

    def get_mesh(self):
        return self.mesh, self.mesh_shift, self.time_symmetry, self.mesh_symmetry

    def set_is_dos_mode(self, is_dos_mode):
        self.is_dos_mode = is_dos_mode

    def get_is_dos_mode(self):
        return self.is_dos_mode

    def set_dos_range(self, dos_min, dos_max, dos_step):
        self.dos_range = { 'min':  dos_min,
                           'max':  dos_max }
        self.omega_step = dos_step

    def get_dos_range(self):
        dos_range = { 'min': self.dos_range['min'],
                      'max': self.dos_range['max'],
                      'step': self.omega_step }
        return dos_range

    def set_pdos_indices(self, indices):
        self.pdos_indices = indices

    def get_pdos_indices(self):
        return self.pdos_indices

    def set_is_thermal_properties(self, is_thermal_properties):
        self.is_thermal_properties = is_thermal_properties

    def get_is_thermal_properties(self):
        return self.is_thermal_properties

    def set_thermal_property_range(self, tmin, tmax, tstep):
        self.tmax = tmax
        self.tmin = tmin
        self.tstep = tstep

    def get_thermal_property_range(self):
        return { 'min':  self.tmin,
                 'max':  self.tmax,
                 'step': self.tstep }

    def set_is_thermal_displacements(self, is_thermal_displacements):
        self.is_thermal_displacements = is_thermal_displacements

    def get_is_thermal_displacements(self):
        return self.is_thermal_displacements

    def set_is_thermal_distances(self, is_thermal_distances):
        self.is_thermal_distances = is_thermal_distances

    def get_is_thermal_distances(self):
        return self.is_thermal_distances

    def set_thermal_atom_pairs(self, atom_pairs):
        self.thermal_atom_pairs = atom_pairs

    def get_thermal_atom_pairs(self):
        return self.thermal_atom_pairs

    def set_anime_band_index(self, band_index):
        self.anime_band_index = band_index

    def get_anime_band_index(self):
        return self.anime_band_index

    def set_anime_amplitude(self, amplitude):
        self.anime_amplitude = amplitude

    def get_anime_amplitude(self):
        return self.anime_amplitude

    def set_anime_division(self, division):
        self.anime_division = division

    def get_anime_division(self):
        return self.anime_division

    def set_anime_qpoint(self, qpoint):
        self.anime_qpoint = qpoint

    def get_anime_qpoint(self):
        return self.anime_qpoint

    def set_anime_shift(self, shift):
        self.anime_shift = shift
    
    def get_anime_shift(self):
        return self.anime_shift

    def set_anime_type(self, anime_type):
        self.anime_type = anime_type
    
    def get_anime_type(self):
        return self.anime_type

    def set_modulation(self, modulation):
        self.modulation = modulation

    def get_modulation(self):
        return self.modulation

        
class PhonopyConfParser( ConfParser ):
    def __init__( self, filename=None, options=None, option_list=None ):
        ConfParser.__init__( self, filename, options, option_list )
        self.__read_options()
        self.__parse_conf()
        self.settings = PhonopySettings()
        self.__set_settings()

    def __read_options( self ):
        for opt in self.option_list:
            if opt.dest=='is_dos_mode':
                if self.options.is_dos_mode:
                    self.confs['dos'] = '.true.'

            if opt.dest=='is_thermal_properties':
                if self.options.is_thermal_properties:
                    self.confs['tprop'] = '.true.'

            if opt.dest=='is_thermal_displacements':
                if self.options.is_thermal_displacements:
                    self.confs['tdisp'] = '.true.'
                    
            if opt.dest=='is_nomeshsym':
                if self.options.is_nomeshsym:
                    self.confs['mesh_symmetry'] = '.false.'

            if opt.dest=='is_read_force_constants':
                if self.options.is_read_force_constants:
                    self.confs['force_constants'] = 'read'
    
            # Overwrite
            if opt.dest=='is_check_symmetry':
                if self.options.is_check_symmetry: # Dummy 'dim' setting for sym-check
                    self.confs['dim'] = '1 1 1'

    def __parse_conf(self):
        confs = self.confs

        for tag in confs.keys():
            if tag == 'band_points':
                self.set_parameter('band_points',
                                   int(confs['band_points']))

            if tag == 'band':
                bands = []
                for section in confs['band'].split(','):
                    points = [ fracval(x) for x in section.split()]
                    if len(points) % 3 != 0:
                        self.setting_error("BAND is incorrectly set.")
                        break
                    bands.append( np.array(points).reshape(-1, 3) )
                self.set_parameter('band', bands )

            if tag == 'force_constants':
                self.set_parameter('force_constants',
                                     confs['force_constants'])

            if tag == 'qpoints':
                self.set_parameter('qpoints', confs['qpoints'])

            if tag == 'mp':
                vals = [ int(x) for x in confs['mp'].split() ]
                if len(vals) < 3:
                    self.setting_error("Mesh numbers are incorrectly set.")
                self.set_parameter('mesh_numbers', vals[:3])

            if tag == 'mp_shift':
                vals = [ fracval(x) for x in confs['mp_shift'].split()]
                if len(vals) < 3:
                    self.setting_error("MP_SHIFT is incorrectly set.")
                self.set_parameter('mp_shift', vals[:3])
                
            if tag == 'mp_reduce':
                self.set_parameter('mp_reduce', confs['mp_reduce'])

            if tag == 'mesh_symmetry':
                self.set_parameter('mesh_symmetry', confs['mesh_symmetry'])

            # Animation
            if tag == 'anime':
                vals = []
                data = confs['anime'].split()
                if len(data) < 3:
                    self.setting_error("ANIME is incorrectly set.")
                
                if len(data) > 2:
                    vals.append( data[0] ) # band index
                    vals.append( data[1] ) # amplitude
                    vals.append( data[2] ) # number of pictures
                    if len(data) > 5:
                        vals.append( data[3] ) # 
                        vals.append( data[4] ) # shift
                        vals.append( data[5] ) #
                self.set_parameter('anime', vals)

            if tag == 'anime_type':
                if ( confs['anime_type'] == 'arc' or
                     confs['anime_type'] == 'v_sim' or
                     confs['anime_type'] == 'poscar' or
                     confs['anime_type'] == 'xyz' or
                     confs['anime_type'] == 'jmol' ):
                    self.set_parameter('anime_type', confs['anime_type'])
                else:
                    self.setting_error("%s is not available for ANIME_TYPE tag." % tags['anime_type'])

            # Modulation
            if tag == 'modulation':
                modulation = {}
                modulation['q'] = [ 0, 0, 0 ]
                modulation['dimension'] = [ 1, 1, 1 ]
                modulation['argument'] = 90
                mod_list = confs['modulation'].split(',')
                is_header_ok = True

                header = mod_list[0].split() 
                if len( header ) == 6 or len( header ) == 7:
                    q = [ fracval( x ) for x in header[:3] ]
                    dimension = [ int( x ) for x in header[3:6] ]
                    modulation['q'] = q
                    modulation['dimension'] = dimension
                    if len( header ) == 7:
                        modulation['argument'] = float( header[6] )
                else:
                    self.setting_error("MODULATION tag is wrongly set.")
                    is_header_ok = False

                if is_header_ok and len( mod_list ) > 1:
                    vals = []
                    for phonon_mode in mod_list[1:]:
                        pair = [ x for x in phonon_mode.split() ]
                        if len( pair ) > 2 or len( pair ) == 0:
                            self.setting_error("MODULATION tag is wrongly set.")
                        elif len( pair ) == 1:
                            vals.append([ int( pair[0] ), 1.0 ])
                        else:
                            vals.append([ int( pair[0] ), float( pair[1] ) ])

                    modulation['modulations'] = vals
                    self.set_parameter('modulation', modulation)

            # DOS
            if tag == 'pdos':
                vals = []
                for sum_set in confs['pdos'].split(','):
                    indices = [ int(x) for x in sum_set.split() ]
                    vals.append(indices)
                self.set_parameter('pdos', vals)

            if tag == 'dos':
                self.set_parameter('dos', confs['dos'])

            if tag == 'dos_range':
                vals = [ float(x) for x in confs['dos_range'].split() ]
                self.set_parameter('dos_range', vals)

            # Thermal properties
            if tag == 'tprop':
                self.set_parameter('tprop', confs['tprop'])

            # Thermal displacement
            if tag == 'tdisp':
                self.set_parameter('tdisp', confs['tdisp'])

            # Thermal distance
            if tag == 'tdistance':
                atom_pairs = []
                for atoms in confs['tdistance'].split(','):
                    pair = [ int(x)-1 for x in atoms.split() ]
                    if len(pair) == 2:
                        atom_pairs.append( pair )
                    else:
                        self.setting_error("TDISTANCE is incorrectly specified.")
                if len(atom_pairs) > 0:
                    self.set_parameter('tdistance', atom_pairs )
            

    def __set_settings( self ):
        ConfParser.set_settings( self )
        tags = self.tags
    
        # Is force constants written or read ?
        if tags.has_key('force_constants'):
            if tags['force_constants'] == 'write':
                self.settings.set_is_force_constants("write")
            elif tags['force_constants'] == 'read':
                self.settings.set_is_force_constants("read")

        # Mesh
        if tags.has_key('mesh_numbers'):
            self.settings.set_run_mode('mesh')
            if tags.has_key('mp_shift'):
                shift = tags['mp_shift']
            else:
                shift = [0.,0.,0.]
    
            time_symmetry = True
            if tags.has_key('mp_reduce'):
                if tags['mp_reduce'] == '.false.':
                    time_symmetry = False
    
            mesh_symmetry = True
            if tags.has_key('mesh_symmetry'):
                if tags['mesh_symmetry'] == '.false.':
                    mesh_symmetry = False
    
            self.settings.set_mesh(tags['mesh_numbers'], shift, time_symmetry, mesh_symmetry)
    
        # band mode
        if tags.has_key('band'):
            if tags.has_key('band_points'):
                npoints = tags['band_points'] - 1
            else:
                npoints = 50
                
            self.settings.set_run_mode('band')
            bands = []
            
            for band_path in tags['band']:
                nd = len( band_path )
                for i in range( nd - 1):
                    diff = ( band_path[i+1] - band_path[i]) / npoints
                    band = [ band_path[i].copy()]
                    q = np.zeros(3)
                    for j in range(npoints):
                        q += diff
                        band.append( band_path[i] + q)
                    bands.append(band)
            self.settings.set_bands(bands)
    
        # Q-points mode
        if tags.has_key('qpoints'):
            if tags['qpoints'] == '.true.':
                self.settings.set_run_mode('qpoints')
    
        # Anime mode
        if tags.has_key('anime_type'):
            self.settings.set_anime_type( tags['anime_type'] )
    
        if tags.has_key('anime'):
            self.settings.set_run_mode('anime')
            anime_type = self.settings.get_anime_type()
            if anime_type=='v_sim':
                qpoints = [ fracval(x) for x in tags['anime'][0:3] ]
                self.settings.set_anime_qpoint( qpoints )
            else:
                self.settings.set_anime_band_index( int( tags['anime'][0] ) )
                self.settings.set_anime_amplitude( float( tags['anime'][1] ) )
                self.settings.set_anime_division( int( tags['anime'][2] ) )
            if len( tags['anime'] ) == 6:
                self.settings.set_anime_shift( 
                    [ fracval(x) for x in tags['anime'][3:6] ] )
    
        # Modulation mode
        if tags.has_key('modulation'):
            self.settings.set_run_mode('modulation')
            self.settings.set_modulation( tags['modulation'] )
    
        # DOS
        if tags.has_key('dos_range'):
            dos_min =  tags['dos_range'][0]
            dos_max =  tags['dos_range'][1]
            dos_step = tags['dos_range'][2]
            self.settings.set_dos_range(dos_min, dos_max, dos_step)
            self.settings.set_is_dos_mode(True)
    
        if tags.has_key('dos'):
            if tags['dos'] == '.true.':
                self.settings.set_is_dos_mode(True)
    
        if tags.has_key('pdos'):
            self.settings.set_pdos_indices(tags['pdos'])
            self.settings.set_is_eigenvectors(True)
            self.settings.set_is_dos_mode(True)
    
        # Thermal properties
        if tags.has_key('tprop'):
            if tags['tprop'] == '.true.':
                self.settings.set_is_thermal_properties(True)
    
        # Thermal displacement
        if tags.has_key('tdisp'):
            if tags['tdisp'] == '.true.':
                self.settings.set_is_thermal_displacements(True)
                self.settings.set_is_eigenvectors(True)
                self.settings.set_mesh_symmetry(False)
    
        # Thermal distance
        if tags.has_key('tdistance'): 
            self.settings.set_is_thermal_distances(True)
            self.settings.set_is_eigenvectors(True)
            self.settings.set_mesh_symmetry(False)
            self.settings.set_thermal_atom_pairs(tags['tdistance'])
    
