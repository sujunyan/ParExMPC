% This script compares the performance of ALADIN and ADMM in terms of 
% as the table shown in the paper. 
% More exmaples can be found on https://github.com/ferreau/mpcBenchmarking.
% --------------------------------------------------------------

addpath("examples");
addpath("@RRLBMPC");

% caseName = 'helicopter'; 
caseName = 'toyExample'; 
% caseName = 'robotArm'; 
cons_mul = 1;      % constraint multiplier
switch caseName
    case 'toyExample'
        prob = example_toyExample;
        % prob.ni = 30;
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
                    'Cx',prob.Cx, 'dx', prob.dx, 'Cu', prob.Cu, 'du', prob.du, ...
                    'ADMM_sigma', 5e-1, 'N', prob.ni, 'delta', 1e-2, 'rho', 1e-6);
                    % 'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul, 'par_flag', true, 'par_threshold', 20);
    case 'robotArm'
        prob = example_robotArm;
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
                    'Cx',prob.Cx, 'dx', prob.dx, 'Cu', prob.Cu, 'du', prob.du, ...
                    'ADMM_sigma', 5e-1, 'N', prob.ni);
end
mpc0 = mpc0.init; 

x0 = prob.x0;
%mpc0.N = 20;
mpc0.maxiter = 10;
tol = 1e-4;

%% Simulate in closed loop -----------------------------
nsim = 100;
res_dict = {};
nx = mpc0.nx; nu = mpc0.nu; N = mpc0.N;

% 
sim_tol = 1e-4;
% method_vec = ["ALADIN"];
% method_vec = ["ADMM"];
method_vec = ["ADMM", "ALADIN"];
% for method = ["ADMM"] %"ALADIN"]
for method = method_vec
    fprintf("Simulation for method: %s \n", method);
    x_vec = [prob.x0];
    u0_vec = [];
    J_vec = [0];
    
    time = 0; % time used by QP solver
    for i = 1:nsim
        fprintf(" Simulation: %d/%d \n", i, nsim);
        x0 =  x_vec(:,end);
        mpc0 = mpc0.updateX0(x0);
        max_iter = 5;
        if i == 1
            % z1 = zeros(nx * N, 1);
            z1 = kron(ones(N,1), x0);
            z2 = zeros(nu * N, 1);
            lam = zeros(nx * N, 1);
            max_iter = 100;
        end
        tic;
        if method == "ADMM"
            [z1, z2, lam, u0] = mpc0.ADMM_solve(z1, z2, lam, max_iter);
        elseif method == "ALADIN"
            [z1, z2, lam, u0] = mpc0.ALADIN_solve(z1, z2, lam, max_iter);
        end
        elapsed = toc;
        time = time + elapsed;
        xn = mpc0.A*x0 + mpc0.B*u0;
        xn_sim = z1(1:nx);
        xn ;
        z2 ;
        
        
        x_vec = [x_vec, xn];
        u0_vec = [u0_vec, u0];
        Jn = J_vec(end) + x0'*mpc0.Q*x0 + u0'*mpc0.R*u0;
        J_vec = [J_vec, Jn];
        if norm(xn - x0, Inf) < sim_tol
            fprintf(" Converged at step %d \n", i);
            break;
        end
        
    end
    res_dict.(method) = struct("x_vec", x_vec, "J_vec", J_vec, "time", time, "u0_vec", u0_vec, "isim", i);
end

for method = method_vec
    fprintf("Method: %s, Total time: %.2f seconds with isim=%d\n", method, res_dict.(method).time, res_dict.(method).isim);
end

% Report the results

% Plot state trajectories
figure(1);
for method = method_vec
    x_vec = res_dict.(method).x_vec;
    nx = size(x_vec, 1); 
    tspan = 0:(size(x_vec, 2) - 1);
    for kk = 1:nx
        subplot(nx, 1, kk);
        hold on;
        plot(tspan, x_vec(kk, :), 'LineWidth', 1.5, 'DisplayName', method);
        title(sprintf("State x_%d Trajectory", kk));
        xlabel("Time Step");
        ylabel(sprintf("x_%d", kk));
        grid on;
        legend show;
    end
end
hold off;

% Plot running cost
figure(2);
hold on
for method = method_vec
    J_vec = res_dict.(method).J_vec;
    tspan = 0:(length(J_vec) - 1);
    plot(tspan, J_vec, 'LineWidth', 1.5, 'DisplayName', method);
    title("Running Cost J");
    xlabel("Time Step");
    ylabel("Cost");
    grid on;
    legend show;
end
hold off;

% Plot control input
figure(3);
for method = method_vec
    u0_vec = res_dict.(method).u0_vec;
    tspan = 0:(size(u0_vec, 2) - 1);
    for kk = 1:nu
        subplot(nu, 1, kk);
        hold on;
        plot(tspan, u0_vec(kk, :), 'LineWidth', 1.5, 'DisplayName', method);
        title(sprintf("Control Input u_0^%d Trajectory", kk));
        xlabel("Time Step");
        ylabel(sprintf("u_%d", kk));
        grid on;
        legend show;
    end
end
hold off;
