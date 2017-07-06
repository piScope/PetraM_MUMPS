##   Makefile
##
##   default variable setting
##
DIRS=etc lib bin sbin share

MAKE=$(shell which make)
PYTHON=$(shell which python)
PYTHONCONFIG=$(shell which python-config)
CC=$(shell which gcc)
SWIG=$(shell which swig)
SWIGFLAG = -Wall -c++ -python

INSTALLDIR=$(shell echo $(PetraM))
ifeq ($(INSTALLDIR),)
   INSTALLDIR := /usr/local/PetraM
endif
BINDIR=$(INSTALLDIR)/bin
LIBDIR=$(INSTALLDIR)/lib


WHOLE_ARCHIVE = --whole-archive
NO_WHOLE_ARCHIVE = --no-whole-archive

MFEM=/usr/local/mfem-3.3
MFEMLIB = mfem
MFEMINCDIR = $(MFEM)
MFEMLNKDIR = $(MFEM)

MFEMSER=/usr/local/mfem-3.3ser
MFEMSERLIB = mfem
MFEMSERINCDIR = $(MFEMSER)
MFEMSERLNKDIR = $(MFEMSER)

HYPRE=/usr/local/hypre-2.11.0
HYPRELIB = HYPRE

METIS5=/usr/local/
METIS5LIB = metis
METIS5LNKDIR = $(METIS5)/lib/

PORD = /usr/local/lib
PORDLIB = pord
PORDLNKDIR = $(PORD)

SCOTCH = /usr/local/lib
SCOTCHLIB = 
SCOTCHLNKDIR =
#SCOTCHINCDIR =

PARMETIS = /usr/local/
PARMETISLIB = parmetis
PARMETISLNKDIR = $(PARMETIS)/lib
#PARMETISINCLIB = 

#
LAPACK = /usr/local/
LAPACKLNKDIR = $(LAPACK)/lib
LAPACKLIB = lapack

#SCALAPACK (for mumps_solve)
SCALAPACK=/usr/local/lib
SCALAPACKLIB = scalapack

# HYPRE
HYPRE 	 = /usr/local/hypre-2.11.0
HYPREINCDIR = $(HYPRE)/include
HYPRELNKDIR = $(HYPRE)/lib

#BLAS (for mumps_solve)
BLAS=/usr/local/lib
BLASLIB = blas
BLASLIBA = 
BLASDIR = $(BLAS)

#PTHREAD (for mumps_solve)
PTHREAD=/usr/local/lib
PTHREADLIB = pthread

#MPI
MPIINCDIR= /opt/local/include/mpich-mp         #mpi.h
MPICHINCDIR    = /opt/local/include/mpich-mp
MPICHLNKDIR    = /opt/local/lib/mpich-mp
MPILIB = mpi
MPICC = mpicc
MPICXX = mpicxx
MPIFC = mpifort
MPIFL = mpifort
MPI4PYINCDIR = $(shell $(PYTHON) -c "import mpi4py;print mpi4py.get_include()")

#numpy
NUMPYINCDIR = $(shell $(PYTHON) -c "import numpy;print numpy.get_include()")

MUMPS = /usr/local/include
MUMPSLIBDIR = /usr/local/lib
MUMPSINCDIR1 = $(MUMPS)/include
MUMPSINCDIR2 = $(MUMPS)/src

OUTC    = -o 
OPTF    = -O  -DALLOW_NON_INIT
OPTL    = -O 
OPTC    = -O
NOCOMPACTUNWIND = 

# MKL
MKL = 
include ./Makefile.local

MFEMINC  = -I$(MFEMINCDIR)
MFEMSERINC  = -I$(MFEMSERINCDIR)
HYPREINC = -I$(HYPREINCDIR)
HYPRELNK = -L$(HYPRELNKDIR) -l$(HYPRELIB)
MPIINC  = -I$(MPIINCDIR)
MPI4PYINC  = -I$(MPI4PYINCDIR)
MUMPSINC = -I$(MUMPSINCDIR1) -I$(MUMPSINCDIR2)
MUMPSCOMMONLIBA = $(MUMPSLIBDIR)/libmumps_common.a
ZMUMPSLIBA = $(MUMPSLIBDIR)/libzmumps.a
SMUMPSLIBA = $(MUMPSLIBDIR)/libsmumps.a
CMUMPSLIBA = $(MUMPSLIBDIR)/libcmumps.a
DMUMPSLIBA = $(MUMPSLIBDIR)/libdmumps.a

export 

default: cext so
.PHONEY:all install

##
build/setup_local.py: Makefile.local
	mkdir -p build
	$(PYTHON) scripts/write_setup_local.py
cext: build/setup_local.py
	$(MAKE) -C petram/ext cext
cxx: build/setup_local.py
	$(MAKE) -C petram/ext cxx
so:
	$(PYTHON) setup.py build
install:
	mkdir -p $(INSTALLDIR)
	for dir in $(INSTALLDIRS); do mkdir -p $(INSTALLDIR)/$$dir; done
	$(PYTHON) setup.py install --prefix=$(INSTALLDIR)	
clean:
	$(MAKE) -C petram/ext clean
	rm -rf build/*

