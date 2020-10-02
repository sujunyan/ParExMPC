function problem = example_toyExample

problem.A = [0.7115 -0.4345;
             0.4345  0.8853];
 
problem.B = [0.2173;
             0.0573];

problem.umin = -5;
problem.umax = 5;
problem.ymin = -inf;
problem.ymax = inf;

problem.Q  = 10*eye(2);  % Q matrix, symmetric, positive semi-definite
problem.R  = 1;          % R matrix, symmetric, postive  definite (*want to penalize the control in all directions)
problem.ni = 10;         % N horizon lengt
problem.C = eye(2);

[~,problem.P] = dlqr(problem.A,problem.B, ... 
    problem.Q,problem.R);

problem.x0 = [10;0];

problem.ur = {[0]};
problem.yr = {[0;0]};
problem.xNr = [0;0];
problem.info.name = "toyExample";

end