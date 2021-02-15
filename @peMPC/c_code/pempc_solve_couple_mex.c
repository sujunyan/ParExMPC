/*
    matlab interface of pempc_solve_couple subroutine, in case the MPT function is not available
*/
#include "mex.h"
#include "pempc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *del_lam;
    double *yout;
    double *xi, *yin, *x0;
    // TODO: should add more error checking part.
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexfunction:outputMismatch",
          "The number of output should be 2!");
    }
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexfunction:inoutMismatch",
          "The number of output should be 3!");
    }

    xi = mxGetPr(prhs[0]);
    yin = mxGetPr(prhs[1]);
    x0 = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(lam_size, 1, mxREAL); // del_lam
    plhs[1] = mxCreateDoubleMatrix(y_size, 1, mxREAL);   // yout
    del_lam = mxGetPr(plhs[0]);
    yout = mxGetPr(plhs[1]);
    solve_couple(xi,yin,x0,del_lam,yout);
}