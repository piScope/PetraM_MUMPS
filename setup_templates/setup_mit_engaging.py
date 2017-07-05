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

mumps_link_args =  [libporda, mumpscommonliba, dmumpsliba, smumpsliba, cmumpsliba,
                    zmumpsliba]

include_dirs = [hypreincdir, mfemincdir, mpichincdir, numpyincdir,
                mpi4pyincdir, numpyincdir, mumpsincdir1, mumpsincdir2]
library_dirs = [hyprelnkdir, pordlnkdir, metis5lnkdir, parmetislnkdir, scalapacklnkdir,  mpichlnkdir]
libraries    = [hyprelib, pordlib, parmetislib, metis5lib, scalapacklib, mpilib, blaslib, pthreadlib]

library_dirs = [x for x in library_dirs if len(x) != 0]
libraries = [x for x in libraries if len(x) != 0]

ext_modules = []

for kk, name in enumerate(modules):
   if kk == 0:
       extra_link_args = mumps_link_args + [name+'.a']
   else:
       extra_link_args = [name+'.a']

   extra_text = [x for x in
                 ['-Wl', whole_archive] +  extra_link_args + 
                 [no_whole_archive] if x != '']

   extra_link_args =  [','.join(extra_text)]
   extra_link_args =  ['-shared-intel', mkl] + extra_link_args + [nocompactunwind]
   extra_link_args =  [x for x in extra_link_args if len(x) != 0]   

   ext_modules.append(Extension(proxy_names[name],
                        sources=sources[name],
                        extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],   
                        extra_link_args = extra_link_args,
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
