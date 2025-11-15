function problem = example_toyExample

    problem.A = [0.7115 -0.4345;
                 0.4345  0.8853];
    
    problem.B = [0.2173;
                 0.0573];

    problem.umin = -5;
    problem.umax = 5;
    problem.dmin = -inf;
    problem.dmax = inf;

    problem.Q = 10*eye(2);  % Q matrix, symmetric, positive semi-definite
    problem.R = 1;          % R matrix, symmetric, positive  definite (*want to penalize the control in all directions)
    problem.ni = 10;        % N horizon length
    problem.C = eye(2);

    [~,problem.P] = dlqr(problem.A,problem.B, ... 
        problem.Q,problem.R);

    problem.x0 = [10;0];

    problem.ur = {[0]};
    problem.yr = {[0;0]};
    problem.xNr = [0;0];
    problem.info.name = "toyExample";

    % Add for RRLB MPC -------------------
    problem.Cx = [eye(2); -eye(2)];  problem.dx = [11; 11; -11; -11];
    problem.Cu = [eye(1); -eye(1)];  problem.du = [5; -5];

end