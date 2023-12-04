% define the problem ------
A = [0.7115,-0.4345; 0.4345, 0.8853];
B = [0.2173;0.0573];
Q = [10,0;0,10];
R = 1;
umin = -5; umax = 5;
N = 3;
x0 = [10;0];

% define the problem in MPT toolbox ------
model = LTISystem('A',A,'B',B);
model.u.min = umin;
model.u.max = umax;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
model.x.with('terminalPenalty');
model.x.terminalPenalty = model.LQRPenalty;
N = 5;
mpt_ctrl = MPCController(model, N);
mpt_ectrl = mpt_ctrl.toExplicit;

%% Comment/Uncomment required option:

%% OPTION 1: Translate MPT-based LTI system to ParExMPC-based MPC controller:
mpc_ctrl = mpt2parexmpc( model, 'N', N );

%% OPTION 2: Translate MPT-based non-explicit MPC controller to ParExMPC-based MPC controller:
mpc_ctrl = mpt2parexmpc( mpt_ctrl );

%% OPTION 3: Translate MPT-based explicit MPC controller to ParExMPC-based MPC controller:
mpc_ctrl = mpt2parexmpc( mpt_ectrl ) ;

%% OPTION 4: Translate MPT-based non-explicit MPC controller to ParExMPC-based MPC controller:
mpc_ctrl = mpt2parexmpc( mpt_ctrl, 'cons_mul',0.5,'mptSolver','mpqp' );

%% OPTION 5: Translate MPT-based LTI system to ParExMPC-based MPC controller:
mpc_ctrl = mpt2parexmpc( mpt_ctrl, 'N', N, 'cons_mul',0.5,'mptSolver','mpqp' );

%% OPTION 6: Translate MPT-based LTI system to ParExMPC-based MPC controller:
mpc_ctrl = model_mpt2parexmpc( model, N );

%% OPTION 7: Translate MPT-based non-explicit MPC controller to ParExMPC-based MPC controller:
mpc_ctrl = controller_mpt2parexmpc( mpt_ctrl );

% use the ParExMPC interface ------
% call help peMPC to see available optional parameters.
mpc_ctrl = mpc_ctrl.build;

%% OPTION 7: Translate MPT-based non-explicit MPC controller to ParExMPC-based MPC controller:
mpc0 = controller_mpt2parexmpc( mpt_ctrl );

% use the ParExMPC interface ------
% call help peMPC to see available optional parameters.
mpc0 = mpc0.build;

% start of MPC simulation ------
maxiter = 5;
Nsim = 100;
[nx,nu] = size(B);
% the tolerance in the MPC iteration
tol = 1e-4;
X = [x0];
U = [];
for k = 1:Nsim
    % solve the first problem with a large maximum iteration
    if (k == 1)
        lam = zeros((N+1)*nx,1);
        z = zeros(N*(nx+nu) + nx ,1);
        [z,lam,u0] = peMPC_controller_mex(x0,100,tol,z,lam);
    else
        [z,lam,u0] = peMPC_controller_mex(x0,maxiter,tol,z,lam);
    end
    % simulate one step of the system dynamics
    x0 = A*x0+B*u0;
    X = [X,x0];
    U = [U,u0];
end
%% Show figures:
% System states
figure, hold on, box on, xlabel('k'), ylabel('x')
stairs([0:Nsim],X(1,:),'-')
stairs([0:Nsim],X(2,:),'--')
stairs([0,Nsim],[0;0],'k:')
legend('x_{1}','x_{2}','origin')
% Control inputs
figure, hold on, box on, xlabel('k'), ylabel('u')
stairs([0:Nsim-1],U,'-')
stairs([0,Nsim],[umin;umin],'k:')
stairs([0,Nsim],[umax;umax],'k.-')
legend('u','u_{min}','u_{max}')