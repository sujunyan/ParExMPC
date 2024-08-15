% define the problem ------
A = [0.7115,-0.4345; 0.4345, 0.8853];
B = [0.2173;0.0573];
Q = [10,0;0,10];
R = 1;
[~,P] = dlqr(A,B,Q,R);
umin = -5; umax = 5;
N = 3;
x0 = [10;0];

% use the ParExMPC interface ------
% call help peMPC to see available optional parameters.
mpc0 = peMPC(A,B,Q,R,P,'umin',umin,'umax',umax,'N',N, 'par_flag',true, 'par_threshold', 20);
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