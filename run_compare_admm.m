% This script compares the performance of pempc and osqp in terms of 
% as the table shown in the paper. 
% More exmaples can be found on https://github.com/ferreau/mpcBenchmarking.
% --------------------------------------------------------------
% mpt_init;
addpath("examples");
addpath("@RRLBMPC");

% caseName = 'helicopter'; 
caseName = 'toyExample'; 
cons_mul = 1;      % constraint multiplier
ospq_flag = true;  % flag to determine wether user want to use osqp
switch caseName
    case 'toyExample'
        prob = example_toyExample;
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
                    'Cx',prob.Cx, 'dx', prob.dx, 'Cu', prob.Cu, 'du', prob.du);
                    % 'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul, 'par_flag', true, 'par_threshold', 20);
        mpc0 = mpc0.init; 
end

x0 = prob.x0;
%mpc0.N = 20;
mpc0.maxiter = 10;
tol = 1e-4;

%% Simulate in closed loop -----------------------------
nsim = 20;
res_dict = {};
nx = mpc0.nx; nu = mpc0.nu; N = mpc0.N;

for method = ["ADMM", "ALADIN"]
    fprintf("Simulation for method: %s \n", method);
    sim_tol = 1e-4;
    x_vec = [x0]
    J_vec = [0];
    time = 0; % time used by QP solver
    for i = 1:nsim
        fprintf(" Simulation: %d/%d \n", i, nsim);
        x0 =  x_vec(:,end);
        mpc0 = mpc0.updateX0(x0);
        if i == 1
            z1 = zeros(nx * N, 1);
            z2 = zeros(nu * N, 1);
            lam = zeros( nx * N, 1);
        end
        if method == "ADMM"
            tic;
            [z1, z2, lam, u0] = mpc0.ADMM_one_iteration(z1, z2, lam);
            elapsed = toc;
            res_dict.(method).time = res_dict.(method).time + elapsed;
        end
        xn = mpc0.A*x0 + mpc0.B*u0;
        x_vec = [x_vec, xn];
        Jn = J_vec(end) + x0'*mpc0.Q*x0 + u0'*mpc0.R*u0;
        J_vec = [J_vec, Jn];
        
    end
    res_dict.(method) = struct("xL", x_vec, "JL", J_vec, "time", time);
end


