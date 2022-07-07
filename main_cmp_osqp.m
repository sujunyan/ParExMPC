% This script compares the performance of pempc and osqp in terms of 
% as the table shown in the paper. 
% More exmaples can be found on https://github.com/ferreau/mpcBenchmarking.
% --------------------------------------------------------------
% mpt_init;
addpath("examples");
addpath("@peMPC");

caseName = 'toyExample'; 
cons_mul = 1;       % constraint multiplier
ospq_flag = true;  % flag to determine wether user want to use osqp
switch caseName
case 'robotArm'
    prob = example_robotArm;
    cons_mul = 0.85; % constraint multiplier
    mpc0 = peMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,'C',prob.C,'D',prob.D,'ur',prob.ur{1}...
             ,'xr',prob.yr{end},'xNr',prob.xNr,'dmin',prob.dmin,'dmax',prob.dmax ...
                ,'umin',prob.umin,'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul);
case 'toyExample'
    prob = example_toyExample;
    mpc0 = peMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,...
                'umin',prob.umin,'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul);
case 'helicopter'
    prob = example_helicopter;
    cons_mul = 0.95; % constraint multiplier
    mpc0 = peMPC(prob.A,prob.B,prob.Q,prob.R,prob.P,'C',prob.C,'ur',prob.ur{1}...
             ,'xr',prob.yr{end},'xNr',prob.xNr,'dmin',prob.dmin,'dmax',prob.dmax ...
                ,'umin',prob.umin,'umax',prob.umax,'N',prob.ni,'cons_mul',cons_mul);
end

x0 = prob.x0;
%mpc0.N = 20;
mpc0 = mpc0.build;
mpc0.maxiter = 10;
tol = 1e-4;
%% setup problems for osqp 
if ospq_flag
osqp_pro = osqp;
osqp_pro.setup(mpc0.QP_H, mpc0.QP_g, mpc0.QP_A, mpc0.QP_l, mpc0.QP_u, ...
                    'warm_start', true,'verbose',false,'eps_abs',tol,'eps_rel',tol,'max_iter',1e5,'polish',true);
end

%% Simulate in closed loop -----------------------------
nsim = 500;
xqpL = [x0];
xL = [x0];
JqpL = [0];
JL = [0];
sim_tol = 1e-4;
qpTime = 0; % time used by QP solver
peTime = 0; % time used by peMPC
cnt = 0;
for i = 1:nsim
    fprintf("%d-th simulation\n",i);
    x0qp = xqpL(:,end);
    x0 = xL(:,end);
    % update the large matrix for 
    mpc0.QP_l(1:mpc0.nx) = -mpc0.A*x0qp;
    mpc0.QP_u(1:mpc0.nx) = -mpc0.A*x0qp;
    if (i == 1) % solve the first problem with osqp
        if ospq_flag
            osqp_pro.update('l',mpc0.QP_l,'u',mpc0.QP_u);
            osqp_res = osqp_pro.solve();
            u0qp = osqp_res.x(1:mpc0.nu);
        end
        z = zeros(mpc0.N*(mpc0.nx+mpc0.nu) + mpc0.nx,1);
        lam = zeros((mpc0.N+1)*mpc0.nx,1);
        [z,lam,u0] = peMPC_controller_mex(x0,100,tol,z,lam);
    else
        % solve osqp -----------------
        if ospq_flag
            tstart1 = tic;
            osqp_pro.update('l',mpc0.QP_l,'u',mpc0.QP_u);
            osqp_res = osqp_pro.solve();
            qpTime = qpTime + toc(tstart1);
            u0qp = osqp_res.x(1:mpc0.nu);
        end
        % solve pempc ----------------------
        tstart2 = tic;
        [z,lam,u0] = peMPC_controller_mex(x0,mpc0.maxiter,tol,z,lam);
        peTime = peTime + toc(tstart2);
    end
    if ospq_flag
        fprintf("u0 err is %e x0 err is %e\n",norm(u0-u0qp,inf),norm(mpc0.C*x0-mpc0.C*x0qp,inf));
        xnqp = prob.A*x0qp + prob.B*u0qp;
    end
    xn = prob.A*x0+prob.B*u0;
    cons_vio_flag = sum(xn>prob.dmax) + sum(xn<prob.dmin); % check if it violate the state constraint
    if cons_vio_flag
        warning("state constraint violate for pempc, consider set smaller cons_mul");
        break;
    end
    if ospq_flag
        xqpL = [xqpL,xnqp];
        JqpL = [JqpL JqpL(end)+objective(prob,u0qp,x0qp)];
    end
    xL = [xL,xn];
    JL = [JL JL(end)+objective(prob,u0,x0)];
    if (objective(prob,u0,x0) < sim_tol)
        cnt = cnt + 1;
        if (cnt >= 20) % add some robustness
            break;
        end
    else
        cnt = 0;
    end
end
nsim0 = length(JL);
if ospq_flag
    fprintf("QPMPC used %f s, %f s per iter \npeMPC used %f s, %f s per iter\n"...
            ,qpTime,qpTime/nsim0,peTime,peTime/nsim0);
    fprintf("The ratio qpTime/peTime is %f\n",qpTime/peTime);
else
    qpTime = NaN;
    JqpL = [NaN];
end
nx = size(x0,1);
tspan = 0:length(JL)-1;
%% plot -----------------------------------------------
figure(1);
for kk = 1:nx
    subplot(nx,1,kk);
    if ospq_flag
        plot(tspan,xL(kk,:),tspan,xqpL(kk,:));
        legend("pempc","osqp");
    else
        plot(tspan,xL(kk,:));
    end
    title(sprintf("x_%d", kk));
end


figure(2);
if ospq_flag
    plot(tspan,JL,tspan,JqpL);
    legend("pempc","osqp");
else
    plot(tspan,JL);
end
title("Running Cost J");


fprintf("The suboptimality compared to MPC: Jpe/Jqp is %f\n",JL(end)/JqpL(end));
fprintf("%16s\t nx\t nu\t N\t maxIter\t peAvgTime\t qpAvgTime\t speedupRatio\t costRatio\t ",'name');
fprintf("costError\t constraint_mul\t");
fprintf("\n");
fprintf("%16s\t %d\t %d\t %d\t %7d\t %9.3e\t %9.3e\t %12.4f\t %9.4f\t "...
        ,prob.info.name,mpc0.nx,mpc0.nu,mpc0.N,mpc0.maxiter,peTime/nsim0,qpTime/nsim0,qpTime/peTime,JL(end)/JqpL(end));
fprintf("%9.3e\t %9.4f\t",(JL(end)-JqpL(end))/JqpL(end), cons_mul);
fprintf("\n");

function obj = objective(prob,u0,x0)
    ur = prob.ur{end};
    xr = pinv(prob.C)* prob.yr{end};
    obj = (u0-ur)'* prob.R * (u0-ur) + (x0-xr)'*prob.C'*prob.Q*prob.C*(x0-xr);
end