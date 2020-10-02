/*
    matlab interface of pempc_get_control. 
*/
#include "mex.h"
#include "pempc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // TODO: should add more error checking part.
    char err_str[100];
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexfunction:outputMismatch",
          "The number of output should be 3!");
    }
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MATLAB:mexfunction:inoutMismatch",
          "The number of output should be 3!");
    }

    if (mxGetM(prhs[0]) != nx || mxGetN(prhs[0]) != 1){
        sprintf(err_str,"The dimension of x0 should be %dx1",nx);
        mexErrMsgIdAndTxt("MATLAB:mexfunction:dimensionMismatch",
          err_str);
    }

    if (mxGetM(prhs[3]) != y_size || mxGetN(prhs[3]) != 1){
        sprintf(err_str,"The dimension of yin should be %dx1",y_size);
        mexErrMsgIdAndTxt("MATLAB:mexfunction:dimensionMismatch",
          err_str);
    }

    if (mxGetM(prhs[4]) != lam_size || mxGetN(prhs[4]) != 1){
        sprintf(err_str,"The dimension of lamin should be %dx1",lam_size);
        mexErrMsgIdAndTxt("MATLAB:mexfunction:dimensionMismatch",
          err_str);
    }

    double *x0 = mxGetPr(prhs[0]);
    double *max_iter_ptr = mxGetPr(prhs[1]);
    size_t max_iter = *max_iter_ptr; // convert from double to int
    double *tol = mxGetPr(prhs[2]);
    double *yin = mxGetPr(prhs[3]);
    double *lam_in = mxGetPr(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(y_size, 1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(lam_size, 1, mxREAL);   
    plhs[2] = mxCreateDoubleMatrix(nu, 1, mxREAL);   
    double *yout = mxGetPr(plhs[0]);
    double *lam_out = mxGetPr(plhs[1]);
    double *u0 = mxGetPr(plhs[2]);
    pempc_get_control(x0,max_iter,*tol,yin,lam_in,
                        yout,lam_out,u0);
}