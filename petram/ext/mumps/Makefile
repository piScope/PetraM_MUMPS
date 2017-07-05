#
#  This file is part of MUMPS 5.0.1, released
#  on Thu Jul 23 17:08:29 UTC 2015
#
#MUMPSDIR = /Users/shiraiwa/src/MUMPS_5.0.1
#MUMPSINC =  -I. -I$(MUMPSDIR)/include -I$(MUMPSDIR)/src
#libdir = /opt/local/lib

#OUTC    = -o 
#OPTF    = -O  -DALLOW_NON_INIT
#OPTL    = -O 
#OPTC    = -O
#CC = mpicc
#CXX = mpicxx
#FC = mpifort
#FL = mpifort

#MFEMINC  = -I/Users/shiraiwa/src/mfem-3.1/
#MPI4PYINC   = -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/mpi4py/include
#HYPREDIR = /Users/shiraiwa/src/hypre-2.11.0/src/hypre
#HYPREINC = -I$(HYPREDIR)/include

IFILE = hypre_to_mumps.i mumps_solve.i
ALLCXX = $(IFILE:.i=._wrap.cxx)

.PHONEY:allso cxx clean

default: mumps_solve.a  hypre_to_mumps.a
all: mumps_solve.a  hypre_to_mumps.a
cext: mumps_solve.a  hypre_to_mumps.a 

cxx: $(ALLCXX)
%.o: %.cpp
	$(MPICXX) -fPIC $(OPTC) $(INCS) $(HYPREINC) $(MFEMINC) $(MPIINC) $(MUMPSINC) -c $^ $(OUTC)$@
%.a: %.o
	ar cruv $@ $^
	ranlib $@
%._wrap.cxx: %.i
	$(SWIG) $(SWIGFLAG) $(MFEMINC) $(MPI4PYINC) $(MUMPSINC) $<
cleancxx:
	rm *.cxx
clean:
	rm -f *.o  *.a *.so
