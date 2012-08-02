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
from phonopy.structure.cells import get_angles, get_cell_parameters, get_cell_matrix
from phonopy.structure.atoms import Atoms
from phonopy.interface.vasp import write_vasp
from phonopy.units import VaspToTHz

class Animation:
    def __init__( self,
                  qpoint,
                  dynamical_matrix,
                  primitive,
                  shift=None ):

        dynamical_matrix.set_dynamical_matrix( qpoint )
        self.eigenvalues, self.eigenvectors = \
            np.linalg.eigh( dynamical_matrix.get_dynamical_matrix() )
        self.qpoint = qpoint
        self.positions = primitive.get_scaled_positions()
        self.symbols = primitive.get_chemical_symbols()
        self.masses = primitive.get_masses()
        self.lattice = primitive.get_cell()
        if not shift==None:
            self.positions = ( self.positions + shift ) % 1
            
    def __set_cell_oriented( self ):
        # Re-oriented lattice xx, yx, yy, zx, zy, zz
        self.angles = get_angles( self.lattice )
        self.cell_params = get_cell_parameters( self.lattice )
        a, b, c = self.cell_params
        alpha, beta, gamma = self.angles
        self.lattice_oriented = get_cell_matrix( a, b, c, alpha, beta, gamma ) 
        self.positions_oriented = \
            self.__set_orientation( np.dot( self.positions, self.lattice ) )

    # For the orientation, see get_cell_matrix
    def __set_orientation( self, vec_cartesian ):
        return np.dot( np.dot( vec_cartesian, np.linalg.inv( self.lattice ) ),
                       self.lattice_oriented )

    def __set_displacements( self, band_index ):
        u = []
        for i, e in enumerate( self.eigenvectors[ :, band_index ] ):
            u.append( e / np.sqrt( self.masses[ i//3 ] ) )

        self.displacements = np.array( u ).reshape( -1, 3 )

    def write_v_sim( self, factor=VaspToTHz, filename="anime.ascii" ):
        self.__set_cell_oriented()
        lat = self.lattice_oriented
        q = self.qpoint
        text  = "# Phonopy generated file for v_sim 3.6\n"
        text += "%15.9f%15.9f%15.9f\n" % ( lat[0,0], lat[1,0], lat[1,1] )
        text += "%15.9f%15.9f%15.9f\n" % ( lat[2,0], lat[2,1], lat[2,2] )
        for s, p in zip( self.symbols, self.positions_oriented ):
            text += "%15.9f%15.9f%15.9f %2s\n" % ( p[0], p[1], p[2], s )

        for i, val in enumerate( self.eigenvalues ):
            if val > 0:
                omega = np.sqrt( val )
            else:
                omega = -np.sqrt( -val )
            self.__set_displacements( i )
            text += "#metaData: qpt=[%f;%f;%f;%f \\\n" % (
                q[0], q[1], q[2], omega * factor )
            for u in self.__set_orientation( self.displacements ):
                text += "#; %f; %f; %f; %f; %f; %f \\\n" % (
                    u[0].real, u[1].real, u[2].real,
                    u[0].imag, u[1].imag, u[2].imag )
            text += "# ]\n"
        w = open(filename, 'w')
        w.write( text )
        w.close()

    def write_arc( self, band_index, amplitude=1, num_div=20,
                   filename="anime.arc" ):
        self.__set_cell_oriented()
        self.__set_displacements( band_index - 1 )
        displacements = self.__set_orientation( self.displacements )

        a, b, c = self.cell_params
        alpha, beta, gamma = self.angles

        text = ""
        text += "!BIOSYM archive 3\n"
        text += "PBC=ON\n"
    
        for i in range( num_div ):
            text += "                                                                        0.000000\n"
            text += "!DATE\n"
            text += "%-4s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n" % (
                "PBC", a, b, c, alpha, beta, gamma)
            positions = self.positions_oriented \
                + ( displacements * np.exp( 2j * np.pi / num_div * i ) ).imag * amplitude
            for j, p in enumerate( positions ):
                text += "%-5s%15.9f%15.9f%15.9f CORE" % (
                    self.symbols[j], p[0], p[1], p[2] )
                text += "%5s%3s%3s%9.4f%5s\n" % (
                    j+1, self.symbols[j], self.symbols[j], 0.0, j+1 )
                
            text += "end\n"
            text += "end\n"
        
        w = open(filename, 'w')
        w.write( text )
        w.close()
            
    def write_xyz_jmol( self,
                        amplitude=10,
                        factor=VaspToTHz,
                        filename="anime.xyz_jmol" ):
        self.__set_cell_oriented()
        text = ""
        for i, val in enumerate( self.eigenvalues ):
            if val > 0:
                freq = np.sqrt( val )
            else:
                freq = -np.sqrt( -val )
            self.__set_displacements( i )
            displacements = self.__set_orientation( self.displacements ) * amplitude
            text += "%d\n" % len(self.symbols)
            text += "q %s , b %d , f %f " % ( str(self.qpoint), i+1, freq * factor )
            text += "(generated by Phonopy)\n" 
            for s, p, u in zip(
                self.symbols, self.positions_oriented, displacements ):
                text += "%-3s  %22.15f %22.15f %22.15f  " % ( s, p[0], p[1], p[2] )
                text += "%15.9f %15.9f %15.9f\n" % ( u[0].real, u[1].real, u[2].real)
        w = open(filename, 'w')
        w.write( text )
        w.close()

    def write_xyz( self, band_index, amplitude=1, num_div=20,
                   factor=VaspToTHz, filename="anime.xyz" ):
        self.__set_cell_oriented()
        freq = self.eigenvalues[band_index - 1]
        self.__set_displacements( band_index - 1 )
        displacements = self.__set_orientation( self.displacements )
        text = ""
        for i in range( num_div ):
            text += "%d\n" % len(self.symbols)
            text += "q %s , b %d , f %f , " % (
                str(self.qpoint), band_index, freq * factor )
            text += "div %d / %d " % ( i, num_div )
            text += "(generated by Phonopy)\n"
            positions = self.positions_oriented \
                + ( displacements * np.exp( 2j * np.pi / num_div * i ) ).imag * amplitude
            for j, p in enumerate( positions ):
                text += "%-3s %22.15f %22.15f %22.15f\n" % (
                    self.symbols[j], p[0], p[1], p[2] )
        w = open(filename, 'w')
        w.write( text )
        w.close()


    def write_POSCAR( self, band_index, amplitude=1, num_div=20,
                      filename="APOSCAR" ):
        self.__set_displacements( band_index - 1 )
        for i in range( num_div ):
            positions = np.dot( self.positions, self.lattice ) \
                + ( self.displacements * np.exp( 2j * np.pi / num_div * i ) ).imag * amplitude
            atoms = Atoms( cell=self.lattice,
                           positions=positions,
                           masses=self.masses,
                           symbols=self.symbols,
                           pbc=True )
            write_vasp( ( filename+"-%03d" ) % i, atoms, direct=True )


