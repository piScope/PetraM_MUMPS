#include <stdio.h>
#include <iostream>
#include <string>
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

#define ICNTL(I)  icntl[(I)-1]  /* macro s.t. indices match documentation */
#define CNTL(I)   cntl[(I)-1]   /* macro s.t. indices match documentation */
#define INFO(I)   info[(I)-1]   /* macro s.t. indices match documentation */
#define INFOG(I)  infog[(I)-1]  /* macro s.t. indices match documentation */
#define RINFO(I)  rinfo[(I)-1]  /* macro s.t. indices match documentation */
#define RINFOG(I) rinfog[(I)-1] /* macro s.t. indices match documentation */
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
  
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}

  // NZ
  void set_nz32(MUMPS_INT nz){id->nz = nz;}
  void set_nz32_loc(MUMPS_INT nz){id->nz_loc = nz;}
  void set_nz(int64_t nnz){id->nnz = nnz;}    
  void set_nz_loc(int64_t nnz){id->nnz_loc = nnz;}    
  
  void set_a(DMUMPS_REAL *a){id->a = a;}
  void set_a_loc(DMUMPS_REAL *a){id->a_loc = a;}  
  void set_rhs(DMUMPS_REAL *rhs){id->rhs = rhs;}
  void set_nrhs(MUMPS_INT nrhs){id->nrhs = nrhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  // distributed RHS
  void set_irhs_loc(MUMPS_INT lrhs_loc, MUMPS_INT *irhs_loc)
  {
     id->irhs_loc = irhs_loc;
     id->lrhs_loc = lrhs_loc;
  }
  void set_nrhs_lrhs_irhs_rhs_loc(MUMPS_INT nloc_rhs,
                                  MUMPS_INT lrhs_loc,
                                  MUMPS_INT *irhs_loc,
				  DMUMPS_REAL *rhs_loc){
    id->nloc_rhs = nloc_rhs;
    id->lrhs_loc = lrhs_loc;
    id->irhs_loc = irhs_loc;
    id->rhs_loc = rhs_loc;    
  }
  // distributed sol
  void set_sol_loc(DMUMPS_REAL *sol_loc,
		   MUMPS_INT lsol_loc,
		   MUMPS_INT *isol_loc){
       id->sol_loc=sol_loc;
       id->lsol_loc=lsol_loc;
       id->isol_loc=isol_loc;       
  }
  void set_saveparam(const char *prefix, const char *dir){
       strcpy(id->save_prefix, prefix);
       strcpy(id->save_dir, dir);
  }
  void set_oocparam(const char *prefix, const char *dir){
       strcpy(id->ooc_prefix, prefix);
       strcpy(id->ooc_tmpdir, dir);
  }
  DMUMPS_REAL *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  DMUMPS_STRUC_C * get_struct(){return id;};
  int get_info(int i){return id->INFO(i);}
  int get_infog(int i){return id->INFOG(i);}
  DMUMPS_REAL get_rinfo(int i){return id->RINFO(i);}
  DMUMPS_REAL get_rinfog(int i){return id->RINFOG(i);}
  std::string get_version_number(void){return std::string(id->version_number);}        
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
  
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}

  // NZ
  void set_nz32(MUMPS_INT nz){id->nz = nz;}
  void set_nz32_loc(MUMPS_INT nz){id->nz_loc = nz;}
  void set_nz(int64_t nnz){id->nnz = nnz;}    
  void set_nz_loc(int64_t nnz){id->nnz_loc = nnz;}    
  
  void set_a(ZMUMPS_COMPLEX *a){id->a = a;}
  void set_a_loc(ZMUMPS_COMPLEX *a){id->a_loc = a;}  
  void set_rhs(ZMUMPS_COMPLEX *rhs){id->rhs = rhs;}
  void set_nrhs(MUMPS_INT nrhs){id->nrhs = nrhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  // distributed RHS
  void set_irhs_loc(MUMPS_INT lrhs_loc, MUMPS_INT *irhs_loc)
  {
     id->irhs_loc = irhs_loc;
     id->lrhs_loc = lrhs_loc;
  }
  void set_nrhs_lrhs_irhs_rhs_loc(MUMPS_INT nloc_rhs,
                                  MUMPS_INT lrhs_loc,
                                  MUMPS_INT *irhs_loc,
				  ZMUMPS_COMPLEX *rhs_loc){
    id->nloc_rhs = nloc_rhs;
    id->lrhs_loc = lrhs_loc;
    id->irhs_loc = irhs_loc;
    id->rhs_loc = rhs_loc;
  }
  // distributed sol
  void set_sol_loc(ZMUMPS_COMPLEX *sol_loc,
		   MUMPS_INT lsol_loc,
		   MUMPS_INT *isol_loc){
       id->sol_loc=sol_loc;
       id->lsol_loc=lsol_loc;
       id->isol_loc=isol_loc;       
  }
  void set_saveparam(const char *prefix, const char *dir){
       strcpy(id->save_prefix, prefix);
       strcpy(id->save_dir, dir);
  }
  void set_oocparam(const char *prefix, const char *dir){
       strcpy(id->ooc_prefix, prefix);
       strcpy(id->ooc_tmpdir, dir);
  }
  ZMUMPS_COMPLEX *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  ZMUMPS_STRUC_C * get_struct(){return id;};
  int get_info(int i){return id->INFO(i);}
  int get_infog(int i){return id->INFOG(i);}
  ZMUMPS_REAL get_rinfo(int i){return id->RINFO(i);}
  ZMUMPS_REAL get_rinfog(int i){return id->RINFOG(i);}
  std::string get_version_number(void){return std::string(id->version_number);}      
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
  void set_cntl(int i,  SMUMPS_REAL v){id->CNTL(i) = v;}  
  void set_icntl(int i, int v){id->ICNTL(i) = v;}
  SMUMPS_REAL get_cntl(int i){return id->CNTL(i);}  
  int get_icntl(int i){return id->ICNTL(i);}
  void set_job(MUMPS_INT n){id->job = n;}
  
  void set_n(MUMPS_INT n){id->n = n;}
  
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}

  // NZ
  void set_nz32(MUMPS_INT nz){id->nz = nz;}
  void set_nz32_loc(MUMPS_INT nz){id->nz_loc = nz;}
  void set_nz(int64_t nnz){id->nnz = nnz;}    
  void set_nz_loc(int64_t nnz){id->nnz_loc = nnz;}    
  
  void set_a(SMUMPS_REAL *a){id->a = a;}
  void set_a_loc(SMUMPS_REAL *a){id->a_loc = a;}  
  void set_rhs(SMUMPS_REAL *rhs){id->rhs = rhs;}
  void set_nrhs(MUMPS_INT nrhs){id->nrhs = nrhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  // distributed RHS
  void set_irhs_loc(MUMPS_INT lrhs_loc, MUMPS_INT *irhs_loc)
  {
     id->irhs_loc = irhs_loc;
     id->lrhs_loc = lrhs_loc;
  }
  void set_nrhs_lrhs_irhs_rhs_loc(MUMPS_INT nloc_rhs,
                                  MUMPS_INT lrhs_loc,
                                  MUMPS_INT *irhs_loc,
				  SMUMPS_REAL *rhs_loc){
    id->nloc_rhs = nloc_rhs;
    id->lrhs_loc = lrhs_loc;
    id->irhs_loc = irhs_loc;
    id->rhs_loc = rhs_loc;    
  }
  // distributed sol
  void set_sol_loc(SMUMPS_REAL *sol_loc,
		   MUMPS_INT lsol_loc,
		   MUMPS_INT *isol_loc){
       id->sol_loc=sol_loc;
       id->lsol_loc=lsol_loc;
       id->isol_loc=isol_loc;       
  }
  void set_saveparam(const char *prefix, const char *dir){
       strcpy(id->save_prefix, prefix);
       strcpy(id->save_dir, dir);
  }
  void set_oocparam(const char *prefix, const char *dir){
       strcpy(id->ooc_prefix, prefix);
       strcpy(id->ooc_tmpdir, dir);
  }
  SMUMPS_REAL *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  SMUMPS_STRUC_C * get_struct(){return id;};
  int get_info(int i){return id->INFO(i);}
  int get_infog(int i){return id->INFOG(i);}
  SMUMPS_REAL get_rinfo(int i){return id->RINFO(i);}
  SMUMPS_REAL get_rinfog(int i){return id->RINFOG(i);}
  std::string get_version_number(void){return std::string(id->version_number);}    
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
  void set_cntl(int i,  CMUMPS_REAL v){id->CNTL(i) = v;}    
  void set_icntl(int i, int v){id->ICNTL(i) = v;}
  CMUMPS_REAL get_cntl(int i){return id->CNTL(i);}      
  int get_icntl(int i){return id->ICNTL(i);}
  void set_job(MUMPS_INT n){id->job = n;}
  
  void set_n(MUMPS_INT n){id->n = n;}
  
  void set_irn(MUMPS_INT *irn){id->irn = irn;}
  void set_jcn(MUMPS_INT *jcn){id->jcn = jcn;}
  void set_irn_loc(MUMPS_INT *irn){id->irn_loc = irn;}
  void set_jcn_loc(MUMPS_INT *jcn){id->jcn_loc = jcn;}

  // NZ
  void set_nz32(MUMPS_INT nz){id->nz = nz;}
  void set_nz32_loc(MUMPS_INT nz){id->nz_loc = nz;}
  void set_nz(int64_t nnz){id->nnz = nnz;}
  void set_nz_loc(int64_t nnz){id->nnz_loc = nnz;}    
  
  void set_a(CMUMPS_COMPLEX *a){id->a = a;}
  void set_a_loc(CMUMPS_COMPLEX *a){id->a_loc = a;}  
  void set_rhs(CMUMPS_COMPLEX *rhs){id->rhs = rhs;}
  void set_nrhs(MUMPS_INT nrhs){id->nrhs = nrhs;}
  void set_lrhs_nrhs(MUMPS_INT lrhs, MUMPS_INT nrhs){
       id->lrhs = lrhs; 
       id->nrhs = nrhs;
  }
  // distributed RHS
  void set_irhs_loc(MUMPS_INT lrhs_loc, MUMPS_INT *irhs_loc)
  {
     id->irhs_loc = irhs_loc;
     id->lrhs_loc = lrhs_loc;
  }
  void set_nrhs_lrhs_irhs_rhs_loc(MUMPS_INT nloc_rhs,
                                  MUMPS_INT lrhs_loc,
                                  MUMPS_INT *irhs_loc,
				  CMUMPS_COMPLEX *rhs_loc){
    id->nloc_rhs = nloc_rhs;
    id->lrhs_loc = lrhs_loc;
    id->irhs_loc = irhs_loc;
    id->rhs_loc = rhs_loc;    
  }
  // distributed sol
  void set_sol_loc(CMUMPS_COMPLEX *sol_loc,
		   MUMPS_INT lsol_loc,
		   MUMPS_INT *isol_loc){
       id->sol_loc=sol_loc;
       id->lsol_loc=lsol_loc;
       id->isol_loc=isol_loc;       
  }
  void set_saveparam(const char *prefix, const char *dir){
       strcpy(id->save_prefix, prefix);
       strcpy(id->save_dir, dir);
  }
  void set_oocparam(const char *prefix, const char *dir){
       strcpy(id->ooc_prefix, prefix);
       strcpy(id->ooc_tmpdir, dir);
  }
  CMUMPS_COMPLEX *get_rhs(void){return id->rhs;}    
  int set_ictrl(int i){return id->ICNTL(i);}
  CMUMPS_STRUC_C * get_struct(){return id;};
  int get_info(int i){return id->INFO(i);}
  int get_infog(int i){return id->INFOG(i);}
  CMUMPS_REAL get_rinfo(int i){return id->RINFO(i);}
  CMUMPS_REAL get_rinfog(int i){return id->RINFOG(i);}
  std::string get_version_number(void){return std::string(id->version_number);}  
};
  
} /* end of namespace */

int libmumps_solve_example_dist(MPI_Comm comm);
int libmumps_solve_example(MPI_Comm comm);

