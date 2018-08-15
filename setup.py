"""
setup.py file for SWIG example
"""
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

print('building mumps solver interface')

## 
setup_dir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
mumps_solve_incdir = os.path.join(setup_dir, "mumps_solve")
mumps_solve_dir = ''

from distutils.core import *
from distutils      import sysconfig


modules= ["mumps_solve", ]

sdir = "petram/ext/mumps/"
sources = {name: [sdir + name + "_wrap.cxx"] for name in modules}

proxy_names = {name: '_'+name for name in modules}

#mumps_link_args =  [libporda, mumpscommonliba, dmumpsliba, smumpsliba, cmumpsliba,
#                    zmumpsliba]


#include_dirs = [mumps_solve_incdir, mpichincdir, numpyincdir,
#                mpi4pyincdir, mumpsincdir, mumpssrcdir]
import numpy
numpyincdir = numpy.get_include()

import mpi4py
mpi4pyincdir = mpi4py.get_include()

mumps_inc_dir = os.getenv("MUMPS_INC_DIR")
mpi_inc_dir = os.getenv("MPI_INC_DIR")

include_dirs = [mumps_solve_incdir, numpyincdir, mpi4pyincdir,
                mumps_inc_dir, mpi_inc_dir]

#lib_list = ["pord", "parmetis", "metis5", "scalapack",  "blas"]
lib_list = []
library_dirs = [os.getenv("MUMPS_SOLVE_DIR")]
#libraries = ["smumps", "dmumps", "cmumps", "zmumps", "mumps_common"]
libraries = ["mumps_solve"]
for lib in lib_list:
    if eval(lib) != "":
        print lib, eval(lib)
        library_dirs.append(eval(lib+ 'lnkdir'))
        libraries.append(eval(lib+'lib'))
        
mkl = os.getenv("MKL")

ext_modules = []
for kk, name in enumerate(modules):
   extra_link_args = []
   #extra_link_args = [sdir + name+'.a']
   '''
   if kk == 0:
       extra_link_args = mumps_link_args + [sdir + name+'.a']
   else:
       extra_link_args = [sdir + name+'.a']

   if whole_archive != '':
       extra_link_args = ['-Wl', whole_archive] +  extra_link_args + [no_whole_archive]
       extra_link_args =  [x for x in extra_link_args if len(x) != 0]   
       extra_link_args =  [','.join(extra_text)]
   '''
   if mkl != '':
       extra_link_args =  ['-shared-intel', mkl] + extra_link_args
   #if nocompactunwind != '':        
   #    extra_link_args.extend([nocompactunwind])
   #extra_link_args =  ['-fopenmp']+[x for x in extra_link_args if len(x) != 0]
   #extra_link_args =  [x for x in extra_link_args if len(x) != 0]   


   ext_modules.append(Extension(proxy_names[name],
                        sources=sources[name],
                        extra_compile_args = ['-DSWIG_TYPE_TABLE=PyMFEM'],   
                        extra_link_args = extra_link_args,
                        include_dirs = include_dirs,
                        library_dirs = library_dirs,
                        libraries = libraries ))

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
