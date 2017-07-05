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

mumps_extra_link_args =  ["-Wl,-whole-archive",
#                    "/home/shiraiwa/lib/libscalapack.a",
                    "/home/shiraiwa/src/MUMPS_5.0.1/lib/libpord.a",
                    mumpscommonliba, dmumpsliba, 
                    smumpsliba, cmumpsliba, zmumpsliba,
                    '/home/shiraiwa/src/scotch_6.0.4/lib/libptesmumps.a',
                    '/home/shiraiwa/src/scotch_6.0.4/lib/libptscotch.a',
                    '/home/shiraiwa/src/scotch_6.0.4/lib/libptscotcherr.a',
                    '/home/shiraiwa/src/scotch_6.0.4/lib/libscotch.a', "-Wl,-no-whole-archive",]
extra_link_args = []

include_dirs = [hypreincdir, mfemincdir, mpichincdir, numpyincdir,
                mpi4pyincdir, numpyincdir, mumpsincdir1, mumpsincdir2]
#library_dirs = [mfemlnkdir, hyprelnkdir, pordlnkdir, parmetislnkdir, mpichlnkdir]
#libraries    = [mfemlib, hyprelib, pordlib, parmetislib, scalapacklib, mpilib, blaslib, pthreadlib]
library_dirs = [hyprelnkdir, pordlnkdir, metis5lnkdir, parmetislnkdir, mpichlnkdir]
libraries    = [hyprelib, pordlib, parmetislib, metis5lib, scalapacklib, lapacklib, mpilib, blaslib, pthreadlib, 'z']


ext_modules = []
for name in modules:
    if name == 'mumps_solve':
         extra_link_args0 = mumps_extra_link_args + [name+'.a']             
    else:
         extra_link_args0 = extra_link_args + [name+'.a']
    ext_modules.append(Extension(proxy_names[name],
                        sources=sources[name],
                        extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],   
#                        extra_link_args = extra_link_args + [name+'.a', "-Wl,-no_compact_unwind"],
                        extra_link_args = extra_link_args0,
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        libraries = libraries ))



setup (name = 'mumps_c_example',
       version = '0.1',
       author      = "S.Shiraiwa",
       description = """C_EXAMPLE wrapper""",
       ext_modules = ext_modules,
       py_modules = modules, 
       )
