MUMPS standalone

This folder includes sample program to test MUMPS and MUMPS wrapper in
standalone enviroment for debuggin.


n
BUILD STEP)
   MUMPS source is assume to exist at the same level as PetraM_MUMPS

   $ export MPICC=mpicc
   $ export MPIFC=mpif90
   $ export MPIFL=mpif90
   $ export OMPFCFLAG=-fopenmp
   $ export OMPLINKFLAG=-fopenmp

   $ make

z_solvetest
   cd dataset_z
   ../z_solvetest



