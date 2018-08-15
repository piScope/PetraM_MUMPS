#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mumps_c_types.h"
#include "dmumps_c.h"
#include "zmumps_c.h"
#include "smumps_c.h"
#include "cmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define JOB_ANALYSIS 1 
#define JOB_FACTORIZATION 2
#define JOB_SOLVE 3
#define JOB_1_2 4
#define JOB_2_3 5
#define JOB_1_2_3 6
#define USE_COMM_WORLD -987654

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I) cntl[(I)-1] /* macro s.t. indices match documentation */
/* No outputs */
/*  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;*/
/* Call the MUMPS package. */


namespace libmumps_solve
{
//   dmumps
class DMUMPS
{
 private:
  DMUMPS_STRUC_C *id;
  bool  initialized;
 public:
  DMUMPS(int par=1, int sym =0, MPI_Comm comm_c = MPI_COMM_NULL);
   ~DMUMPS();
  void run();
  void finish();
  void set_cntl(int i,  DMUMPS_REAL v){id->CNTL(i) = v;}  
  void set_icntl(int i, int v){id->ICNTL(i) = v;}
  DMUMPS_REAL get_cntl(int i){return id->CNTL(i);}  
  int get_icntl(int i){return id->ICNTL(i);}
  void set_job(MUMPS_INT n){id->job = n;}  
  void set_n(MUMPS_INT n){id->n = n;}
  void set_nz(MUMPS_INT nz){id->nz = nz;}
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_nz_loc(MUMPS_INT nz){id->nz_loc = nz;}  
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}  
  void set_a(DMUMPS_REAL *a){id->a = a;}
  void set_a_loc(DMUMPS_REAL *a){id->a_loc = a;}  
  void set_rhs(DMUMPS_REAL *rhs){id->rhs = rhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  DMUMPS_REAL *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  DMUMPS_STRUC_C * get_struct(){return id;};  
};

//   zmumps  
class ZMUMPS
{
 private:
  ZMUMPS_STRUC_C *id;
  bool  initialized;  
 public:
  ZMUMPS(int par=1, int sym = 0,  MPI_Comm comm_c = MPI_COMM_NULL);
   ~ZMUMPS();
  void run();
  void finish();
  void set_cntl(int i,  ZMUMPS_REAL v){id->CNTL(i) = v;}
  void set_icntl(int i, int v){id->ICNTL(i) = v;}
  ZMUMPS_REAL get_cntl(int i){return id->CNTL(i);}
  int get_icntl(int i){return id->ICNTL(i);}
  void set_job(MUMPS_INT n){id->job = n;}  
  void set_n(MUMPS_INT n){id->n = n;}
  void set_nz(MUMPS_INT nz){id->nz = nz;}
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_nz_loc(MUMPS_INT nz){id->nz_loc = nz;}  
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}  
  void set_a(ZMUMPS_COMPLEX *a){id->a = a;}
  void set_a_loc(ZMUMPS_COMPLEX *a){id->a_loc = a;}  
  void set_rhs(ZMUMPS_COMPLEX *rhs){id->rhs = rhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  ZMUMPS_COMPLEX *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  ZMUMPS_STRUC_C * get_struct(){return id;};  
};

//   smumps  
class SMUMPS
{
 private:
  SMUMPS_STRUC_C *id;
  bool  initialized;  
 public:
  SMUMPS(int par=1, int sym = 0, MPI_Comm comm_c = MPI_COMM_NULL);
   ~SMUMPS();
  void run();
  void finish();
  void set_icntl(int i, int v){id->ICNTL(i) = v;}
  int get_icntl(int i){return id->ICNTL(i);}
  void set_job(MUMPS_INT n){id->job = n;}  
  void set_n(MUMPS_INT n){id->n = n;}
  void set_nz(MUMPS_INT nz){id->nz = nz;}
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_nz_loc(MUMPS_INT nz){id->nz_loc = nz;}  
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}  
  void set_a(SMUMPS_REAL *a){id->a = a;}
  void set_a_loc(SMUMPS_REAL *a){id->a_loc = a;}  
  void set_rhs(SMUMPS_REAL *rhs){id->rhs = rhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  SMUMPS_REAL *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  SMUMPS_STRUC_C * get_struct(){return id;};  
};

//   cmumps  
class CMUMPS
{
 private:
  CMUMPS_STRUC_C *id;
  bool  initialized;  
 public:
  CMUMPS(int par=1, int sym = 0, MPI_Comm comm_c = MPI_COMM_NULL);
   ~CMUMPS();
  void run();
  void finish();
  void set_icntl(int i, int v){id->ICNTL(i) = v;}
  int get_icntl(int i){return id->ICNTL(i);}
  void set_job(MUMPS_INT n){id->job = n;}  
  void set_n(MUMPS_INT n){id->n = n;}
  void set_nz(MUMPS_INT nz){id->nz = nz;}
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_nz_loc(MUMPS_INT nz){id->nz_loc = nz;}  
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}  
  void set_a(CMUMPS_COMPLEX *a){id->a = a;}
  void set_a_loc(CMUMPS_COMPLEX *a){id->a_loc = a;}  
  void set_rhs(CMUMPS_COMPLEX *rhs){id->rhs = rhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  CMUMPS_COMPLEX *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  CMUMPS_STRUC_C * get_struct(){return id;};  
};
  
} /* end of namespace */

int libmumps_solve_example_dist(MPI_Comm comm);
int libmumps_solve_example(MPI_Comm comm);

