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
        prob.ni = 2;
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
                    'Cx',prob.Cx, 'dx', prob.dx, 'Cu', prob.Cu, 'du', prob.du, ...
                    'ADMM_sigma', 5e-1, 'N', prob.ni);
                    % 'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul, 'par_flag', true, 'par_threshold', 20);
        mpc0 = mpc0.init; 
end

x0 = prob.x0;
%mpc0.N = 20;
mpc0.maxiter = 10;
tol = 1e-4;

%% Simulate in closed loop -----------------------------
nsim = 50;
res_dict = {};
nx = mpc0.nx; nu = mpc0.nu; N = mpc0.N;

for method = ["ADMM"] %"ALADIN"]
    fprintf("Simulation for method: %s \n", method);
    sim_tol = 1e-4;
    x_vec = [x0];
    J_vec = [0];
    
    time = 0; % time used by QP solver
    for i = 1:nsim
        fprintf(" Simulation: %d/%d \n", i, nsim);
        x0 =  x_vec(:,end);
        mpc0 = mpc0.updateX0(x0);
        max_iter = 5;
        if i == 1
            z1 = zeros(nx * N, 1);
            z2 = zeros(nu * N, 1);
            lam = zeros( nx * N, 1);
            max_iter = 100;
        end
        if method == "ADMM"
            tic;
            [z1, z2, lam, u0] = mpc0.ADMM_solve(z1, z2, lam, max_iter);
            elapsed = toc;
            time = time + elapsed;
        end
        xn = mpc0.A*x0 + mpc0.B*u0;
        xn_sim = z1(1:nx);
        xn ;
        
        
        x_vec = [x_vec, xn];
        Jn = J_vec(end) + x0'*mpc0.Q*x0 + u0'*mpc0.R*u0;
        J_vec = [J_vec, Jn];
        
    end
    res_dict.(method) = struct("x_vec", x_vec, "J_vec", J_vec, "time", time);
end

% Extract the number of states and simulation time steps
nx = size(x_vec, 1);
tspan = 0:(size(x_vec, 2) - 1);

% Plot state trajectories
figure(1);
for kk = 1:nx
    subplot(nx, 1, kk);
    plot(tspan, x_vec(kk, :), 'LineWidth', 1.5);
    title(sprintf("State x_%d Trajectory", kk));
    xlabel("Time Step");
    ylabel(sprintf("x_%d", kk));
    grid on;
end

% Plot running cost
figure(2);
plot(tspan, J_vec, 'LineWidth', 1.5);
title("Running Cost J");
xlabel("Time Step");
ylabel("Cost");
grid on;
