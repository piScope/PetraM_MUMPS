#
#  
#
IFILE = mumps_solve.i
ALLCXX = $(IFILE:.i=._wrap.cxx)

.PHONEY:allso cxx clean

default: mumps_solve.a
all: mumps_solve.a 
cext: mumps_solve.a

cxx: $(ALLCXX)
%.o: %.cpp
	$(MPICXX) -fPIC $(OPTC) $(INCS) $(MPIINC) $(MUMPSINC) -c $^ $(OUTC)$@
%.a: %.o
	ar cruv $@ $^
	ranlib $@
%._wrap.cxx: %.i
	$(SWIG) $(SWIGFLAG) $(MFEMINC) $(MPI4PYINC) $(MUMPSINC) $<
cleancxx:
	rm *.cxx
clean:
	rm -f *.o  *.a *.so
