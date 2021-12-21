function problem = example_robotArm

problem.Ts = 0.004; % [s]

Ac = [  0.00e+0,  1.00e+0,  0.00e+0,  0.00e+0; ...
        0.00e+0, -1.72e+1,  0.00e+0,  0.00e+0; ...
        0.00e+0,  0.00e+0,  0.00e+0,  1.00e+0; ...
        0.00e+0,  1.00e+0,  0.00e+0, -1.61e+1  ];

Bc = [  0.00e+0,  0.00e+0; ...
        2.62e+0,  0.00e+0; ...
        0.00e+0,  0.00e+0; ...
        0.00e+0,  2.48e+0  ];
    
Cc = eye(4);

sys = c2d( ss(Ac,Bc,Cc,[]),problem.Ts,'zoh' );
[A,B,C,D] = ssdata(sys);


problem.A = A;        % compulsory
problem.B = B;        % compulsory
problem.C = C;        % optional (default, identity matrix)
problem.D = D;        % optional (default, zero matrix)

problem.umax =  [ 100; 25 ];            % all constraint fields are optional
problem.umin = -[ 100; 25 ];            % (default values -inf/inf)

problem.dmax =  [ Inf; 1; Inf; 1 ];     % all constraint fields are optional
problem.dmin = -[ Inf; 1; Inf; 1 ];     % (default values -inf/inf)

problem.Q  = diag([ 1.14e+4, 2.24e+1, 1.14e+4, 2.94e+1 ]);  % 
problem.R  = diag([ 2.20e-1, 2.37e-1 ]);                    % 
problem.P = problem.Q;
         
problem.x0 = [ -1; 0; 1; 0 ];
problem.ni = 20;

problem.ur = {[0;0]};
problem.yr = {[0;0;0;0]};
problem.xNr = [0;0;0;0];
problem.info.name = "robotArm";

end