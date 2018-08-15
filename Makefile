##   Makefile
##
##   default variable setting
##
MAKE=$(shell which make)
PYTHON ?= $(shell which python)
PYTHONCONFIG ?= $(shell which python-config)
CC ?= $(shell which gcc)
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

#
#  these are to absorb the difference between linux and macOS
#
# no-compact-unwind (this option is not used anymore...)
NOCOMPACTUNWIND ?=
ifeq ($(UNAME), Darwin)
   NOCOMPACTUNWIND = -Wl,-no_compact_unwind
endif

MPI_INC_DIR ?= $(INSTALLDIR)/include
MUMPS_SOLVE_DIR ?= $(INSTALLDIR)/lib
MUMPS_INC_DIR   ?= $(INSTALLDIR)/include
MKL ?=
CC=${MPICC}
CXX=${MPICXX}

export 

default: so
.PHONEY:all install

##
cxx: 
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

