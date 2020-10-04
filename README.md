# ParExMPC Toolbox 

ParExMPC Toolbox a MATLAB base toolkit for linear MPC with state constraints. In particular, we solve the problem of the following form:
$$
\text{minimize}\; (x_N - x_{r})^T P (x_N - x_{r}) + 
    \sum_{k=0}^{N-1} (x_k-x_r)^T Q (x_k-x_r) + (u_k - u_r)^T R (u_k - u_r) \\
\text{subject to }  x_{k+1} = Ax_k+Bx_k \\
\qquad y_{min} \leq Cx+Du \leq y_{max} \\
u_{min} \leq u\leq u_{max}
$$

## Getting Started

Please follow the instruction on the website and install the [Multi-Parametric Toolbox 3](https://www.mpt3.org/).
One can construct a pempc object in MATLAB via
``` matlab
mpc0 = peMPC(A,B,Q,R,P,{,optional inputs});
mpc0.init;
```
where A,B are system dynamics and Q,R,P are cost parameters in linear MPC setting. To generate c codes, run:
``` matlab
mpc0.getMPTfunc;
mpc0.toC;
mpc0.compile_c;
```
where the generated c code can be found in the folder @peMPC/c_code. Then one can use the mex function for MPC iteration
``` matlab
[z,lam,u0] = pempc_get_control_mex(x0,mpc0.maxiter,tol,z,lam);
```
For a quick start, see the example script `main_cmp_osqp.m`. 

## C code on Embedded system

To see how one can directly use generated c codes on Embedded system, see our Arduino example in the folder examples/arduino_example.

## Reference

Y. Jiang, J. Oravec, B. Houska and M. Kvasnica, "Parallel MPC for Linear Systems with Input Constraints," in IEEE Transactions on Automatic Control, doi: 10.1109/TAC.2020.3020827.
