#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

print('building mumps solver interface')

## first load variables from PyMFEM_ROOT/setup_local.py
root =  os.path.abspath(os.path.join(os.path.realpath(__file__),
                                     '..', 'build'))
print root
sys.path.insert(0, root)
from  setup_local import *

from distutils.core import *
from distutils      import sysconfig

modules= ["mumps_solve", "hypre_to_mumps"]

sdir = "petram/ext/mumps/"
sources = {name: [sdir + name + "_wrap.cxx"] for name in modules}

proxy_names = {name: '_'+name for name in modules}

extra_link_args =  [mumpscommonliba, dmumpsliba, smumpsliba, cmumpsliba,
                    zmumpsliba]
include_dirs = [hypreincdir, mfemincdir, mpichincdir, numpyincdir,
                mpi4pyincdir, numpyincdir, mumpsincdir1, mumpsincdir2]
#library_dirs = [mfemlnkdir, hyprelnkdir, pordlnkdir, parmetislnkdir, mpichlnkdir]
#libraries    = [mfemlib, hyprelib, pordlib, parmetislib, scalapacklib, mpilib, blaslib, pthreadlib]
library_dirs = [hyprelnkdir, pordlnkdir, metis5lnkdir, parmetislnkdir, mpichlnkdir]
libraries    = [hyprelib, pordlib, parmetislib, metis5lib, scalapacklib, mpilib, blaslib, pthreadlib]


ext_modules = [Extension(proxy_names[name],
                        sources=sources[name],
                        extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],   
                        extra_link_args = extra_link_args + [sdir + name+'.a', "-Wl,-no_compact_unwind"],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        libraries = libraries )
              for name in modules]


setup (name = 'PetraM_MUMPS',
       url='https://github.com/piScope/PetraM',       
       version = '0.1',
       description = 'PetraM MUMPS interface', 
       long_description=long_description,       
       author      = "S. Shiraiwa",
       author_email='shiraiwa@psfc.mit.edu',
       license='GNUv3',
       
       classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 2.7',
       ],

       keywords='MFEM physics',
       packages=find_packages(),
       ext_modules = ext_modules,
       py_modules = modules,
)
