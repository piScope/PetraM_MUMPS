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

z_solvetest/c_solvetest: load matrix and rhs from file and solve
   * matrix is distirubted
   * rhs is host

   cd dataset_complex
   ../z_solvetest
   ../c_solvetest   



