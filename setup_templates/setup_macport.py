#!/usr/bin/env python

"""
setup.py file for SWIG example
"""
print('building mumps solver interface')

## first load variables from PyMFEM_ROOT/setup_local.py
import sys
import os
root =  os.path.abspath(os.path.join(os.path.realpath(__file__),
                                     '..', '..', '..', '..'))
sys.path.insert(0, root)
from  setup_local import *

from distutils.core import *
from distutils      import sysconfig

modules= ["mumps_solve", "hypre_to_mumps"]

sources = {name: [name + "_wrap.cxx"] for name in modules}

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
                        extra_link_args = extra_link_args + [name+'.a', "-Wl,-no_compact_unwind"],
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        libraries = libraries )
              for name in modules]


setup (name = 'mumps_c_example',
       version = '0.1',
       author      = "S.Shiraiwa",
       description = """C_EXAMPLE wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
