/*
    matlab interface of pempc_solve_couple subroutine, in case the MPT function is not available
*/
#include "mex.h"
#include "pempc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *xiout;
    double *yin, *lam_in, *x0;
    // TODO: should add more error checking part.
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexfunction:outputMismatch",
          "The number of output should be 1!");
    }
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexfunction:inoutMismatch",
          "The number of output should be 3!");
    }
    yin = mxGetPr(prhs[0]);
    lam_in = mxGetPr(prhs[1]);
    x0 = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(y_size, 1, mxREAL);   // xiout
    xiout = mxGetPr(plhs[0]);
    solve_decouple(yin,lam_in,x0,xiout);
}