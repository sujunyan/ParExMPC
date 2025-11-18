% This script compares the performance of ALADIN and ADMM in terms of 
% as the table shown in the paper. 
% More exmaples can be found on https://github.com/ferreau/mpcBenchmarking.
% --------------------------------------------------------------

addpath("examples");
addpath("@RRLBMPC");

% caseName = 'helicopter'; 
caseName = 'toyExample'; 
caseName = 'robotArm'; 
cons_mul = 1;      % constraint multiplier
prob = feval(['example_', caseName]);

mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
            'Cx',prob.Cx, 'dx', prob.dx, 'Cu', prob.Cu, 'du', prob.du, ...
            'ADMM_sigma', 5e-1, 'N', prob.ni, 'delta', 1e-2, 'rho', 1e-6);
mpc0 = mpc0.init; 
x0 = prob.x0;
mpc0.maxiter = 10;
tol = 1e-4;
nx = mpc0.nx; nu = mpc0.nu; N = mpc0.N;

%% Simulate in closed loop -----------------------------
nsim = 20;
res_dict = {};

% 
sim_tol = 1e-4;
% method_vec = ["ALADIN"];
% method_vec = ["fmincon", "ALADINiter5", "ADMMiter5",  "ADMMiter20", "ADMMiter50" ];
method_vec = ["fmincon", ];
% method_vec = ["fmincon"];
% method_vec = ["ADMM", "ALADIN"];
% for method = ["ADMM"] %"ALADIN"]
for method_str = method_vec
    if method_str ~= "fmincon"
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

% Report the results

% Plot state trajectories
figure(1);
lw = 4.0;
for method = method_vec
    x_vec = res_dict.(method).x_vec;
    nx = size(x_vec, 1); 
    tspan = 0:(size(x_vec, 2) - 1);
    lw = lw - 0.4;
    for kk = 1:nx
        subplot(nx, 1, kk);
        hold on;
        plot(tspan, x_vec(kk, :), 'LineWidth', lw, 'DisplayName', method);
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
lw = 4.0;
for method = method_vec
    lw = lw - 0.4;
    J_vec = res_dict.(method).J_vec;
    tspan = 0:(length(J_vec) - 1);
    plot(tspan, J_vec, 'LineWidth', lw, 'DisplayName', method);
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
