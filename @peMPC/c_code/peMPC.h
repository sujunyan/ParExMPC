
#ifndef __PEMPC_H__
#define __PEMPC_H__

#include "predefined.h"
#include <stddef.h>

#define y_size (N*(nx+nu)+nx)
#define lam_size ((N+1)*nx)

// self-defined linear algebra routine -----------------------
void pempc_aAxpy(const char trans, const double a, const double* A, 
               const double* x, const size_t M, const size_t N, double* y);

void pempc_aAx(const char trans, const double a, const double* A, 
               const double* x, const size_t M, const size_t N, double* y);

void pempc_axpby(const double a, const double *x, const double b, 
                const double *y, const size_t N, double * out);

// main ALADIN algorithm -------------------------
void solve_couple(const double *xi, const double *yin, const double *x0,    
                  double *del_lam, double *yout);
void solve_decouple(const double *yin, const double *lam_in, const double *x0,
                    double *xiout);
void pempc_get_control(const double *x0, size_t max_iter, double tol, 
            const double *yin ,const double *lam_in, 
            double *yout, double *lam_out, double *u0);

// MPT related functions --------------------------------
unsigned long mpt_eval(const double *X, const size_t i_start, double *U, 
const size_t domain, const size_t range, const size_t mpt_nr, const int * mpt_nc, 
const double *A, const double *B, const double *F, const double *G, 
const double *HTB, const double *GTB, const double *FTB);

unsigned long MPT_func_0( double *X, double *U);
unsigned long MPT_func_k( double *X, double *U);
unsigned long MPT_func_N( double *X, double *U);

#endif