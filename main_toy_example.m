mpt_init;
addpath(genpath("@peMPC")); 
addpath("examples");

% define the problem ------
  Q = [10,0;0,10]; 
  R = 1;
  A = [0.7115,-0.4345; 0.4345, 0.8853];
  B = [0.2173;0.0573];
  [~,P] = dlqr(A,B,Q,R);
  umin = -5; umax = 5;
  N = 3;
  x0 = [10;0];
  
  % use the ParExMPC interface ------
  mpc0 = peMPC(A,B,Q,R,P,'umin',umin,'umax',umax,'N',N);
  mpc0 = mpc0.init;                % initialize the pempc solver
  mpc0 = mpc0.getMPTfunc; 	% generate mpt codes
  mpc0.toC;               	% generate c code
  mpc0.compile_c;        
  maxiter = 5;
  % start of MPC simulation ------
  nsim = 100;
  [nx,nu] = size(B);
  tol = 1e-4;               % the tolerance in the MPC iteration
  for i = 1:nsim
    if (i == 1)             % solve the first problem with large maximum iteration
      lam = zeros((N+1)*nx,1);
      z = zeros(N*(nx+nu) + nx ,1);
      [z,lam,u0] = peMPC_controller_mex(x0,100,tol,z,lam);
    else
      [z,lam,u0] = peMPC_controller_mex(x0,maxiter,tol,z,lam);
    end
    x0 = A*x0+B*u0;         % simulate one step of the system dynamics
  end
  