#include <string.h>
#include "pempc.h"
#include <time.h>
#include <float.h>
#include "mpQP_data.h"
#include <assert.h>
#include <omp.h>
//#define USE_BLAS
#ifdef USE_BLAS 
#include "blas.h" // This file is copied from matlab
#endif
// #define USE_OPENMP (PEMPC_N >= 20) change to controlled by the MATLAB API code


// commented out for the use of dualized x0
//static unsigned long MPT_func_0_( double *X, double *U, const size_t istart){
//    return mpt_eval(X,istart,U,nx+nu,nu,zMPT_NR,zMPT_NC,zMPT_A,zMPT_B,zMPT_F,zMPT_G,zMPT_HTB,zMPT_GTB,zMPT_FTB);
//}
static unsigned long MPT_func_N_( double *X, double *U, const size_t istart){
    return mpt_eval(X,istart,U,nx,nx,NMPT_NR,NMPT_NC,NMPT_A,NMPT_B,NMPT_F,NMPT_G,NMPT_HTB,NMPT_GTB,NMPT_FTB);
}
static unsigned long MPT_func_k_( double *X, double *U, const size_t istart){
    return mpt_eval(X,istart,U,nx+nu,nx+nu,kMPT_NR,kMPT_NC,kMPT_A,kMPT_B,kMPT_F,kMPT_G,kMPT_HTB,kMPT_GTB,kMPT_FTB);
}

/* *****************************************************
    some simple operations of linear algebra.
    Some functions are the wrappers of blas. Can be extended to self-implemented vertion in case BLAS library is not availabe.
********************************************************/
// out = a*x + b*y, the vector size 
void pempc_axpby(const double a, const double *x, const double b, 
                const double *y, const size_t N, double * out){
    for (size_t i = 0; i<N; i++){
        out[i] = a*x[i] + b*y[i];
    }
}
/**
 * @brief: calute the result of matrix y = a*A*x
 * @note In case the blas is not available, use self-defined instead, 
   @param: if trans == 'n': y = a*A*x; 
           if trans == 't': y = a*A'*x
 * @param: M,N : the dimension of A is M by N
 * @retval: y: the return array
    TODO: no error-handling currently
*/
void pempc_aAx(const char trans, const double a, const double* A, 
               const double* x, const size_t M, const size_t N, double* y){
    size_t LDA = M;
    static size_t inc = 1;
    static double dzero = 0;
    #ifdef USE_BLAS
    dgemv(&trans,&M,&N,&a,A,&LDA,x,&inc,&dzero,y,&inc);
    #else
    // decide the size of x and y
    size_t nx_,ny_;
    if (trans == 't'){
        nx_ = M; ny_ = N;
    }else if (trans == 'n'){
        nx_ = N; ny_ = M;
    }else{
        assert(0);
    }
    double x_tmp[nx_];
    double y_tmp[nx_];
    memset(y_tmp,0,ny_*sizeof(double));
    memcpy(x_tmp,x,nx_*sizeof(double));
    pempc_aAxpy(trans,a,A,x_tmp,M,N,y_tmp);
    memcpy(y,y_tmp,ny_*sizeof(double));
    #endif // end of USE_BLAS
}

/**
 * @brief: calute the result of matrix y = aAx + y
 * @note In case the blas is not available, use self-defined instead, 
   @param: if trans == 'n': y = a*A*x; 
           if trans == 't': y = a*A'*x
 * @param: M,N : the dimension of A is M by N
 * @retval: y: the return array
    TODO: no error-handling currently
*/
void pempc_aAxpy(const char trans, const double a, const double* A, 
               const double* x, const size_t M, const size_t N, double* y){
    size_t LDA = M;
    static size_t inc = 1;
    static double one = 1;
    #ifdef USE_BLAS
    dgemv(&trans,&M,&N,&a,A,&LDA,x,&inc,&one,y,&inc);
    #else
    // decide the size of x and y
    size_t nx_,ny_;
    double Aij;
    if (trans == 't'){
        nx_ = M; ny_ = N;
    }else if (trans == 'n'){
        nx_ = N; ny_ = M;
    }else{
        assert(0);
    }
    double y_tmp[ny_];
    memcpy(y_tmp,y,ny_*sizeof(double));
    for (size_t iRow = 0; iRow<M; iRow++){
        for (size_t iCol = 0; iCol<N; iCol++){
            Aij = A[iRow + iCol*M];
            if (trans == 'n'){
                y_tmp[iRow] += a* Aij*x[iCol];
            }else{
                y_tmp[iCol] += a* Aij*x[iRow];
            }
        }
    }
    memcpy(y,y_tmp,ny_*sizeof(double));
    #endif // end of USE_BLAS
}

#if 0
// For debug
void printArray(const double *A, size_t N){
    for (size_t i = 0; i<N;i++){
        printf("%f,\t",A[i]);
    }
    printf("\n");
}
#endif

/*****************************************************************
    functions for peMPC controller
*************************************************/

/**
 * @brief  get the MPC control u0 
 * @note   
 * @param  *x0:         nx x 1 vector, the initial state 
 * @param  max_iter:    the max number of iteration
 * @param  tol:         if norm(xi-y) < tol, terminate the process
 * @param  *yin:        (N*(nx+nu)+nx)x1 vector: the initial guess of primal variable y 
 * @param  *lam_in:     (N*nx+nx)x1 vector: the initial guess of dual variable lam 
 * @retval *yout:       (N*(nx+nu)+nx)x1 vector: the shifted primal variable for warm-start in next sampling time
 * @retval *lam_out:    (N*nx+nx)x1 vector: the shifted dual variable for warm-start in next sampling time
 * @retval *u0:         the control input  
 */
void pempc_get_control(const double *x0, size_t max_iter, double tol, 
            const double *yin ,const double *lam_in, 
            double *yout, double *lam_out, double *u0){
    double xi[y_size];
    double y_tmp[y_size];
    double lam_tmp[lam_size];
    double lam_tmp2[lam_size];
    double del_lam[lam_size];
    int ter_flag = 1; // the flag to determine if the tolerence is reached.
    memcpy(y_tmp,yin,sizeof(y_tmp)); // y_tmp = yin;
    memcpy(lam_tmp,lam_in,sizeof(lam_tmp)); // lam_tmp = lam_in;
    for (size_t m = 0; m < max_iter; m++){
        solve_decouple(y_tmp,lam_tmp,x0,xi);
        solve_couple(xi,y_tmp,x0,del_lam,yout);
        pempc_axpby(1,lam_tmp,1,del_lam,lam_size,lam_tmp2); // lam_tmp2 = lam_tmp + del_lam 
        memcpy(lam_tmp,lam_tmp2,sizeof(lam_tmp));
        pempc_axpby(1,xi,-1,yout,y_size,y_tmp);  // y_tmp = xi - yout
        ter_flag = 1;
        for (size_t i = 0; i < y_size; i++){
            if(y_tmp[i] > tol){
                ter_flag = 0;
                break;
            }
        }
        if(ter_flag){break;} // if the termination condition is reached
        memcpy(y_tmp,yout,sizeof(y_tmp)); // y_tmp = yout
    }
    // get the outputs
    memcpy(u0,xi+nx,nu*sizeof(double));
    // shift the primal and dual variables
    memcpy(y_tmp,yout,sizeof(y_tmp)); // y_tmp = yout
    memcpy(yout,y_tmp+nu+nx,(y_size-nx-nu)*sizeof(double)); // yout{1:N-1} = {u1,x2,..,x{N-1},u{N-1},x{N}}
    memset(yout+(N-1)*(nx+nu)+nx,0,(nx+nu)*sizeof(double)); // yout{N} = zeros;
    // dual variable
    memcpy(lam_out,lam_tmp+nx,(lam_size-nx)*sizeof(double));
    memset(lam_out+lam_size-nx,0,nx*sizeof(double)); // lam{N+1} = zeros
}

/**
 * @brief  solve the decouple problems with MPT functions
 * @note   
 * @param  *yin:    (N*(nx+nu)+nx)x1 vector, the primal variable   
 * @param  *lam_in: (N*nx+nx)x1 vector, the dual variable 
 * @param  *x0:     nx x 1 vector, the initial state x0
 * @retval *xiout:  (N*(nx+nu)+nx)x1 vector, the alternating direction variable
 */
void solve_decouple(const double *yin, const double *lam_in, const double *x0,
                    double *xiout){
    double xu_tmp[(N+1)*(nx+nu)];
    double x_tmp[(N+1)*nx];   double x_tmp2[(N+1)*nx];
    double u_tmp[(N+1)*nu];   double u_tmp2[(N+1)*nu];
    const double *yn, *lamn;
    double *xioutn;
    static size_t start_index_array[PEMPC_N];
    const size_t u_array_size = sizeof(double)*nu;
    const size_t x_array_size = sizeof(double)*nx;
    // solve decouple0 --------------------------------------
    // gs0 = -2*R*(ur+y0) + B'*lam0
    #if 0
    pempc_axpby(1,ur,1,yin,nu,u_tmp); // u_tmp = ur+y0
    pempc_aAx('n',-2,R,u_tmp,nu,nu,u_tmp2); // u_tmp2 = -2*R*u_tmp 
    pempc_aAxpy('t',1,B,lam_in,nx,nu,u_tmp2); // u_tmp2 += B'*lam0 ( == gs0)
    memcpy(xu_tmp,x0,x_array_size); 
    memcpy(xu_tmp + nx,u_tmp2,u_array_size);  // input of MPT_func_0 is [x0;gs0]
    MPT_func_0_(xu_tmp,xiout,0);    // xi0 = MPT_func_0(xu_tmp)
    #endif

    // solve decoupleN ------------------------------------
    // gsN = -2*P*(xNr + yN) - lamN
    yn = yin + y_size - nx;
    lamn = lam_in + lam_size - nx;
    xioutn = xiout + y_size - nx;
    pempc_axpby(1,xNr,1,yn,nx,x_tmp); // x_tmp = xNr + yn
    pempc_aAx('n',-2,P,x_tmp,nx,nx,x_tmp2); // x_tmp2 = -2*P*(xNr + yn)
    pempc_axpby(1,x_tmp2,-1,lamn,nx,x_tmp); // x_tmp = x_tmp2 - lamN ( == gsN)
    MPT_func_N_(x_tmp,xioutn,0);  // xiN = MPT_func_N(x_tmp)

    // solve decouplek ---------------------------------
    //clock_t tstart = clock();
    #if USE_OPENMP // to avoid the parallel overhead.
    #pragma omp parallel for private(yn,xioutn,lamn)
    for (size_t n = 0; n < N; n++){
    #else
    for (size_t n = 0; n < N; n++){
    #endif
        double *x_tmpn, *x_tmp2n, *u_tmpn, *u_tmp2n, *xu_tmpn;
        x_tmpn = x_tmp + n*nx;
        x_tmp2n = x_tmp2 + n*nx;
        u_tmpn = u_tmp + n*nu;
        u_tmp2n = u_tmp2 + n*nu;
        xu_tmpn = xu_tmp + n*(nu+nx);
        yn = yin + n*(nx+nu);
        xioutn = xiout + n*(nx+nu);
        lamn = lam_in + n*nx;
        // x part gsn_x = -2*Q*(xr+ynx) - lamn + A'*lam_{n+1}
        pempc_axpby(1,xr,1,yn,nx,x_tmpn); //x_tmpn = xr + yn_x
        pempc_aAx('n',-2,Q,x_tmpn,nx,nx,x_tmp2n); // x_tmp2n = -2*Q*x_tmpn
        pempc_aAxpy('t',1,A,lamn+nx,nx,nx,x_tmp2n); //  x_tmp2n += A'*lam_{n+1}
        pempc_axpby(1,x_tmp2n,-1,lamn,nx,x_tmpn); // x_tmpn = x_tmp2n - lamn (== gsn_x)
        // u part gsn_u = -2*R*(ur+yu) + B'*lam_{n+1}
        pempc_axpby(1,ur,1,yn+nx,nu,u_tmpn); //u_tmpn = ur + yn_u
        pempc_aAx('n',-2,R,u_tmpn,nu,nu,u_tmp2n); // u_tmp2n = -2*R*u_tmpn
        pempc_aAxpy('t',1,B,lamn+nx,nx,nu,u_tmp2n); //  u_tmp2n += B'*lam_{n+1} (== gsn_u)
        // xu_tmpn = gsn
        memcpy(xu_tmpn,x_tmpn,x_array_size);
        memcpy(xu_tmpn+nx,u_tmp2n,u_array_size);
        unsigned long region = MPT_func_k_(xu_tmpn,xioutn,start_index_array[n]); // xioutn = MPT_func_k(xu_tmpn)
        start_index_array[n] = region-1;
        //printf("MPT_func_k region: %ld\n",region);
    }
    //printf("solve decouplek used %f s\n",(double)(clock()-tstart)/CLOCKS_PER_SEC);
}


/**
 * @brief   solve the coupled problem with Riccati based method
 * @note   
 * @param  *xi:     (N*(nx+nu)+nx)x1 vector, the alternating direction variable
 * @param  *yin:    (N*(nx+nu)+nx)x1 vector, the previous primal variable
 * @param  *x0:     nx x 1 vector, the initial state x0
 * @retval  *del_lam:del_lam:(N*nx))x1 vector, the increase of dual variable
 * @retval  *yout:   (N*(nx+nu))x1 vector, the next primal variable
 */
void solve_couple(const double *xi, const double *yin, const double *x0,    
                  double *del_lam, double *yout){
    double p[(N+1)*nx];
    double l[N*nu];
    double qn[nx];  // the buffer array for intermedita usage
    double sn[nu];  // the buffer array for intermedita usage
    double x_tmp[nx]; // buffer of size nx
    double u_tmp[nu]; // buffer of size nu
    double xn[nx]; // 
    double un[nu]; // 
    double ximy[y_size]; // the buffer that store the variable 2*xi - yin
    double *pn, *ximyn, *ln, *Ric_Lam_invn, *Ric_Ln, *Ric_Pn, *yn;  
    double *del_lamn;
    pempc_axpby(2,xi,-1,yin,y_size,ximy); // ximy = 2*xi-yin
    pn = p + nx*N;
    ximyn = ximy + y_size - nx;
    pempc_aAx('n',-2,P,ximyn,nx,nx,pn); // pN = -2*P*ximyN

    // backward sweep
    for (int n = N-1; n>=0; n--){
        ximyn = ximy + n*(nx+nu);
        pempc_aAx('n',-2,Q,ximyn,nx,nx,qn); // qn = -2*Q*tmp
        pempc_aAx('n',-2,R,ximyn+nx,nu,nu,sn); // sn = -2*R*tmp
        ln = l + nu*n; 
        pn = p + nx*n;
        Ric_Lam_invn = Ric_Lam_inv + (nu*nu)*n;
        Ric_Ln = Ric_L + (nx*nu)*n;
        pempc_aAxpy('t',1,B,pn+nx,nx,nu,sn); // sn = sn + B'*p{n+1}
        pempc_aAx('n',1,Ric_Lam_invn,sn,nu,nu,ln); // ln = Ric_Lam_invn * sn
        pempc_aAxpy('t',1,A,pn+nx,nx,nx,qn); // qn = qn + A'*p{n+1}
        pempc_aAxpy('t',-1,Ric_Ln,ln,nu,nx,qn); // qn = qn - Ric_Ln'*ln
        memcpy(pn,qn,sizeof(qn)); // pn = qn
    }
    // forward routine
    memcpy(xn,x0,sizeof(xn)); // xn = x0
    memcpy(del_lam,p,sizeof(xn)); // del_lam{1} = p{1}
    pempc_aAxpy('n',1,Ric_P, xn, nx,nx,del_lam); // del_lam{1} += Ric_P{1} * x0
    for (size_t n =0; n<=N-1;n++){
        ln = l + n*nu;
        Ric_Ln = Ric_L + (nu*nx)*n;
        Ric_Lam_invn = Ric_Lam_inv + (nu*nu)*n;
        del_lamn = del_lam + (n+1)*nx;
        Ric_Pn = Ric_P + n*(nx*nx);
        pn = p + n*nx;
        memcpy(u_tmp,ln,sizeof(u_tmp)); // u_tmp = ln
        pempc_aAxpy('n',1,Ric_Ln,xn,nu,nx,u_tmp); // u_tmp = u_tmp + Ric_Ln * xn
        pempc_aAx('t',-1,Ric_Lam_invn,u_tmp,nu,nu,un);  // un = - Ric_Lam_inv' * u_tmp
        yn = yout + n*(nx+nu);
        memcpy(yn,xn,sizeof(xn));
        memcpy(yn+nx,un,sizeof(un)); // yn = [xn;un]

        pempc_aAx('n',1,A,xn,nx,nx,x_tmp); // x_tmp = A*xn
        pempc_aAxpy('n',1,B,un,nx,nu,x_tmp); // x_tmp += B*un 
        memcpy(xn,x_tmp,sizeof(x_tmp)); // xn = x_tmp;
        memcpy(del_lamn,pn+nx,sizeof(xn)); // del_lam{n+1} = p{n+1}
        pempc_aAxpy('t',1,Ric_Pn + (nx*nx), xn, nx,nx,del_lamn);    // del_lam{n+1}+= Ric_P{n+1}' * xn
    }
    yn = yout + N*(nx+nu);
    memcpy(yn,xn,sizeof(xn)); // yN = xN
}


const double MPT_ABSTOL = 1.000000e-8;
/**
 * @brief mpt evaluation function using sequential search 
 * @note   
 * @retval  *U :   nux1 vector, the optimal solution of the mpqp
 */
unsigned long mpt_eval(const double *X, const size_t i_start, double *U, 
const size_t domain, const size_t range, const size_t mpt_nr, const int * mpt_nc, 
const double *A, const double *B, const double *F, const double *G, 
const double *HTB, const double *GTB, const double *FTB){
    int ix, jx, ic, nc, isinside;
    size_t ireg, abspos, iregmin, region, ireg0;
    double hx, sx, obj, objmin;
    const double *An;
    size_t inc = 1;

    abspos = 0;
    region = 0;
    iregmin = 0;

   /* initialize values of the tie-break function */
    obj = 0;
    objmin = DBL_MAX;

    //printf("start mpt_eval--------------------\n");
    // make the right offset
    for (ireg = 0; ireg < i_start; ireg++){
        //printf("ireg %d,istart is %d\n",ireg,i_start);
        abspos += mpt_nc[ireg];
    }

    //printf("istart is %d\n",i_start);
    for (ireg0=0; ireg0<mpt_nr; ireg0++) {
        ireg = (i_start + ireg0)%mpt_nr;
        //printf("tring ireg %d, NR %d\n",ireg,mpt_nr);
        //ireg = ireg0;
        isinside = 1;
        nc = mpt_nc[ireg];
        for (ic=0; ic<nc; ic++) {
            #if 1
            hx = 0;
            for (ix=0; ix<domain; ix++) {
                hx += A[abspos*domain+ic*domain+ix]*X[ix];
            }
            #else
            An = A+(abspos*domain) + (ic*domain);
            hx = ddot(&domain,An,&inc,X,&inc); 
            #endif
            if ((hx - B[abspos+ic]) > MPT_ABSTOL) {
                /* constraint is violated, continue with next region */
                isinside = 0;
                break;
            } 
        }
        #if 0
        if (isinside==1) {
            /* state belongs to this region, evaluate the tie-breaking function */
            obj = 0;
            for (ix=0; ix<domain; ix++) {
                sx = 0;
                for (jx=0; jx<domain; jx++) {
                    sx += HTB[ireg*domain*domain + ix*domain + jx]*X[jx];
                }
                obj += sx*X[ix];
            }
            for (ix=0; ix<domain; ix++) {
                obj += FTB[ireg*domain + ix]*X[ix];
            }
            obj += GTB[ireg];
            //printf("detect region %d with obj %e min obj %e\n",ireg,obj,objmin);
            if (obj<objmin) {
                objmin = obj;
                region = ireg + 1;
                iregmin = ireg;
            }
        }
        abspos = abspos + mpt_nc[ireg];
        #else
        if (isinside == 1){ // TODO: delete the tiebreak check, might be wrong.
            //printf("detect region %d with isinside %d\n",ireg,isinside);
            //objmin = obj;
            region = ireg + 1;
            iregmin = ireg;
            break;
        }
        if (ireg == mpt_nr -1){
            // reach the end, make the circle back to index 0
            abspos = 0;
        }else{
            abspos += mpt_nc[ireg];
        }
        #endif
    }
    for (ix=0; ix<range; ix++) {
        sx = 0;
        for (jx=0; jx<domain; jx++) {
            sx += F[iregmin*domain*range + ix*domain + jx]*X[jx];
        }
        U[ix] = sx + G[iregmin*range + ix];
    }
    //printf("mpt_eval end-----------------------\n");
    return region;
}

