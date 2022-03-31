#
#  This file is part of MUMPS 5.4.1, released
#  on Tue Aug  3 09:49:43 UTC 2021
#
topdir = ../../MUMPS_5.4.1
mumpsdir = ../../MUMPS_5.4.1
libdir = $(mumpsdir)/lib

default: z

.PHONY: default all s d c z multi clean
.SECONDEXPANSION:

all:	c z s d multi

c:	
z:	z_solvetest
s:	
d:	
multi:

include $(mumpsdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)
LIBSMUMPS = $(libdir)/libsmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)
LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)
LIBZMUMPS = $(libdir)/libzmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)
LIBCMUMPS = $(libdir)/libcmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)


z_solvetest:  $(LIBSMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) z_solvetest.o  $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) -I. -I$(mumpsdir)/include -I$(mumpsdir)/src $(INCS) -c $*.F $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -I$(mumpsdir)/include -I$(mumpsdir)/src $(INCS) -c $*.c $(OUTC)$*.o


clean:
	$(RM) *.o z_solvetest



