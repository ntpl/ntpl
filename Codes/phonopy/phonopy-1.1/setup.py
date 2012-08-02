from distutils.core import setup, Extension
#from setuptools import setup, Extension
import numpy
include_dirs_numpy = [numpy.get_include()]

extension = Extension('phonopy._phonopy',
                      include_dirs = ['c'] + include_dirs_numpy,
                      sources = ['c/_phonopy.c'])

extension_spglib = Extension('phonopy._spglib',
                             include_dirs = ['c/spglib'] + include_dirs_numpy,
                             # extra_compile_args=['-fopenmp'],
                             # extra_link_args=['-lgomp'],
                             sources = ['c/_spglib.c',
                                        'c/spglib/refinement.c',
                                        'c/spglib/cell.c',
                                        'c/spglib/debug.c',
                                        'c/spglib/hall_symbol.c',
                                        'c/spglib/kpoint.c',
                                        'c/spglib/lattice.c',
                                        'c/spglib/mathfunc.c',
                                        'c/spglib/pointgroup.c',
                                        'c/spglib/primitive.c',
                                        'c/spglib/spacegroup.c',
                                        'c/spglib/spg_database.c',
                                        'c/spglib/spglib.c',
                                        'c/spglib/site_symmetry.c',
                                        'c/spglib/sitesym_database.c',
                                        'c/spglib/symmetry.c'] )

setup( name = 'phonopy',
       version = '1.1',
       description = 'This is the phonopy module.',
       author = 'Atsushi Togo',
       author_email = 'atz.togo@gmail.com',
       url = 'http://phonopy.sourceforge.net/',
       packages = ['phonopy',
                   'phonopy.harmonic',
                   'phonopy.interface',
                   'phonopy.hphonopy',
                   'phonopy.cui',
                   'phonopy.phonon',
                   'phonopy.structure'],
       scripts = ['scripts/phonopy',
                  'scripts/phonopy-qha',
                  'scripts/phonopy-FHI-aims',
                  'scripts/bandplot',
                  'scripts/propplot',
                  'scripts/tdplot',
                  'scripts/dispmanager',
                  'scripts/pdosplot'],
       ext_modules = [extension, extension_spglib] )
