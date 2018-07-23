##   Makefile
##
##   default variable setting
##
MAKE=$(shell which make)
PYTHON ?= $(shell which python)
PYTHONCONFIG ?= $(shell which python-config)
CC ?= $(shell which gcc)
SWIG ?= $(shell which swig)
SWIGFLAG = -Wall -c++ -python
UNAME := $(shell uname)

INSTALLDIR=$(shell echo $(PetraM))
ifeq ($(INSTALLDIR),)
   INSTALLDIR := /usr/local/PetraM
endif
BINDIR=$(INSTALLDIR)/bin
LIBDIR=$(INSTALLDIR)/lib

### compilers
MPILIB ?= mpi
MPICC ?= mpicc
MPICXX ?= mpicxx
MPIFC ?= mpifort
MPIFL ?= mpifort

### library locations
USRLOCAL ?= /usr/local

MUMPS ?= $(USRLOCAL)
PORD ?= $(USRLOCAL)

### For the followin libraries,  if variabls are set empty,
### it is not used in build process
METIS5 ?= $(USRLOCAL)
SCOTCH ?= $(USRLOCAL)
PARMETIS ?= $(USRLOCAL)
LAPACK ?= $(USRLOCAL)
SCALAPACK ?= $(USRLOCAL)
BLAS ?= $(USRLOCAL)

#MPI
MPIINCDIR ?= /opt/local/include/mpich-mp         #mpi.h
MPICHINCDIR ?= /opt/local/include/mpich-mp
MPICHLNKDIR ?= /opt/local/lib/mpich-mp


MPI4PYINCDIR = $(shell $(PYTHON) -c "import mpi4py;print mpi4py.get_include()")
NUMPYINCDIR = $(shell $(PYTHON) -c "import numpy;print numpy.get_include()")

#OMP
OMPFLAG ?= ""   # (gcc) "-fopenmp" (intel compiler) "-qopenmp"

OUTC    ?= -o 
OPTF    ?= -O  -DALLOW_NON_INIT
OPTL    ?= -O 
OPTC    ?= -O

# MKL
MKL ?=

#
#  these are to absorb the difference between linux and macOS
#
# no-compact-unwind
NOCOMPACTUNWIND ?=
ifeq ($(UNAME), Darwin)
   NOCOMPACTUNWIND = -Wl,-no_compact_unwind
endif

WHOLE_ARCHIVE = --whole-archive
NO_WHOLE_ARCHIVE = --no-whole-archive
ifeq ($(UNAME), Darwin)
   WHOLE_ARCHIVE = 
   NO_WHOLE_ARCHIVE =
endif

include ./Makefile.local

SCOTCHLIB = 
ifneq ($(SCOTCH), "")
   SCOTCHLIB = scotch
   SCOTCHLNKDIR = $(SCOTCH)/lib
   SCOTCHINCDIR = $(SCOTCH)/include
endif

PARMETISLIB = 
ifneq ($(PARMETIS), "")
   PARMETISLIB = parmetis
   PARMETISLNKDIR = $(PARMETIS)/lib
endif

METIS5LIB =
ifneq ($(METIS5), "")
   METIS5LIB = metis
   METIS5LNKDIR = $(METIS5)/lib/
endif

LAPACKLIB = 
ifneq ($(LAPACK), "")
   LAPACKLIB = lapack
   LAPACKLNKDIR = $(LAPACK)/lib
endif

BLASLIB = 
ifneq ($(BLAS), "")
   BLASLIB = blas
   BLASLNKDIR=$(BLAS)/lib
   BLASINCDIR=$(BLAS)/include
endif

SCALAPACKLIB =
ifneq ($(SCALAPACK), "")
   SCALAPACKLIB = scalapack
   SCALAPACKLNKDIR=$(SCALAPACK)/lib
endif

PORDLIB = pord
PORDLNKDIR = $(PORD)/lib


MUMPSLIBDIR =  $(MUMPS)/lib
MUMPSINCDIR = $(MUMPS)/include
MUMPSSRCDIR = $(MUMPS)/src
MUMPSINC = -I$(MUMPSINCDIR) -I$(MUMPSSRCDIR)
LIBPORDA = $(MUMPSLIBDIR)/libpord.a
MUMPSCOMMONLIBA = $(MUMPSLIBDIR)/libmumps_common.a
ZMUMPSLIBA = $(MUMPSLIBDIR)/libzmumps.a
SMUMPSLIBA = $(MUMPSLIBDIR)/libsmumps.a
CMUMPSLIBA = $(MUMPSLIBDIR)/libcmumps.a
DMUMPSLIBA = $(MUMPSLIBDIR)/libdmumps.a

MPIINC  = -I$(MPIINCDIR)
MPI4PYINC  = -I$(MPI4PYINCDIR)


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

