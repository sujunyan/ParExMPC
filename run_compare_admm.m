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
    case 'robotArm'
        prob = example_robotArm;
        cons_mul = 0.85; % constraint multiplier
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,'C',prob.C,'D',prob.D,'ur',prob.ur{1}...
                 ,'xr',prob.yr{end},'xNr',prob.xNr,'dmin',prob.dmin,'dmax',prob.dmax ...
                    ,'umin',prob.umin,'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul, 'par_flag', true, 'par_threshold', 20);
    case 'toyExample'
        prob = example_toyExample;
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
                    'umin',prob.umin,'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul, 'par_flag', true, 'par_threshold', 20);
    case 'helicopter'
        prob = example_helicopter;
        cons_mul = 0.95; % constraint multiplier
        mpc0 = RRLBMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,'C',prob.C,'ur',prob.ur{1}...
                 ,'xr',prob.yr{end},'xNr',prob.xNr,'dmin',prob.dmin,'dmax',prob.dmax ...
                    ,'umin',prob.umin,'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul, 'par_flag', true, 'par_threshold', 20);
end

x0 = prob.x0;
%mpc0.N = 20;
mpc0 = mpc0.build;
mpc0.maxiter = 10;
tol = 1e-4;

%% Simulate in closed loop -----------------------------
nsim = 500;
xqpL = [x0];
xL = [x0];
JqpL = [0];
JL = [0];
sim_tol = 1e-4;
admmTime = 0; % time used by QP solver
aladinTime = 0; % time used by peMPC
cnt = 0;
for i = 1:nsim
    fprintf(" Simulation: %d/%d \n", i,nsim);
    x0qp = xqpL(:,end);
    x0 = xL(:,end);

    % solve admm ------


    % solve aladin ------
    
end

