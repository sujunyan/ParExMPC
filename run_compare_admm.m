% This script compares the performance of ALADIN and ADMM in terms of 
% as the table shown in the paper. 
% More exmaples can be found on https://github.com/ferreau/mpcBenchmarking.
% --------------------------------------------------------------

addpath("examples");
addpath("@RRLBMPC");

% caseName = 'helicopter'; 
% caseName = 'toyExample'; 
caseName = 'robotArm'; 
cons_mul = 1;      % constraint multiplier
prob = feval(['example_', caseName]);
rho = 1e1;
mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
            'Cx',prob.Cx, 'dx', prob.dx, 'Cu', prob.Cu, 'du', prob.du, ...
            'ADMM_sigma', 5e-1, 'N', prob.ni, 'delta', 1e-2, 'rho', rho);
mpc0 = mpc0.init; 
x0 = prob.x0;
mpc0.maxiter = 10;
tol = 1e-4;
nx = mpc0.nx; nu = mpc0.nu; N = mpc0.N;

%% Simulate in closed loop -----------------------------
nsim = 5;
res_dict = {};

% 
sim_tol = 1e-4;
% method_vec = ["ALADIN"];
% method_vec = ["fmincon", "ALADINiter5", "ADMMiter5",  "ADMMiter20", "ADMMiter50" ];
method_vec = ["quadprog", "fmincon"];
% method_vec = ["fmincon"];
% method_vec = ["ADMM", "ALADIN"];
% for method = ["ADMM"] %"ALADIN"]
for method_str = method_vec
    if (method_str ~= "fmincon") && (method_str ~= "quadprog")
        parts = strsplit(method_str, 'iter');
        method = parts{1};
        max_iter0 = str2num(parts{2});
    else
        method = method_str;
        max_iter0 = 0;
    end
    
    fprintf("Simulation for method: %s \n", method_str);
    res = simulate_one(mpc0, nsim, sim_tol, method, max_iter0, prob.x0); 

    res_dict.(method_str) = res;
end

for method = method_vec
    fprintf("Method: %s, Total time: %.2f seconds with isim=%d\n", method, res_dict.(method).time, res_dict.(method).isim_stop);
end

plot_sim_results(res_dict);