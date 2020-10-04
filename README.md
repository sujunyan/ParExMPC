# ParExMPC Toolbox 

ParExMPC Toolbox a MATLAB base toolkit for linear MPC with state constraints. In particular, we solve the problem of the following form:

<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bequation%7D%0A%5Cbegin%7Barray%7D%7Bll%7D%0A%5Ctext%7Bminimize%7D%5C%3B%26%20%28x_N%20-%20x_%7Br%7D%29%5ET%20P%20%28x_N%20-%20x_%7Br%7D%29%20%2B%20%0A%20%20%20%20%5Csum_%7Bk%3D0%7D%5E%7BN-1%7D%20%28x_k-x_r%29%5ET%20Q%20%28x_k-x_r%29%20%2B%20%28u_k%20-%20u_r%29%5ET%20R%20%28u_k%20-%20u_r%29%20%5C%5C%0A%5Ctext%7Bsubject%20to%7D%5C%3B%26%20%20x_%7Bk%2B1%7D%20%3D%20Ax_k%2BBx_k%20%5C%5C%0A%26y_%7Bmin%7D%20%5Cleq%20Cx%2BDu%20%5Cleq%20y_%7Bmax%7D%20%5C%5C%0A%26u_%7Bmin%7D%20%5Cleq%20u%5Cleq%20u_%7Bmax%7D%0A%5Cend%7Barray%7D%0A%5Cend%7Bequation%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\begin{equation}\begin{array}{ll}\text{minimize}\;& (x_N - x_{r})^T P (x_N - x_{r}) +     \sum_{k=0}^{N-1} (x_k-x_r)^T Q (x_k-x_r) + (u_k - u_r)^T R (u_k - u_r) \\\text{subject to}\;&  x_{k+1} = Ax_k+Bx_k \\&y_{min} \leq Cx+Du \leq y_{max} \\&u_{min} \leq u\leq u_{max}\end{array}\end{equation}" width="714" height="85" />

<!---
$$
\text{minimize}\; (x_N - x_{r})^T P (x_N - x_{r}) + 
    \sum_{k=0}^{N-1} (x_k-x_r)^T Q (x_k-x_r) + (u_k - u_r)^T R (u_k - u_r) \\
\text{subject to }  x_{k+1} = Ax_k+Bx_k \\
\qquad y_{min} \leq Cx+Du \leq y_{max} \\
u_{min} \leq u\leq u_{max}
$$
--->

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
