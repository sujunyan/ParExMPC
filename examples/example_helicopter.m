function [ problem ] = Benchmark_helicopter
% try to load benchmark from mat file
problem.info.name         = 'helicopter';
            
problem.A = [   0.99       0       0.01    0       0       0;...
                0       0.99       0       0.01    0       0;...
                0       0       0.99       0       0       0;...
                0       0       0       0.99       0       0;...
                0.01    0       0       0       0.99       0;...
                0       0.01    0       0       0       0.99];


problem.B = [   0           0;...
                0.0001      -0.0001;...
                0.0019      -0.0019;...
                0.0132      -0.0132;...
                0           0;...
                0           0];
            
problem.umax =  [   3;     3    ];             
problem.umin =  [  -1;    -1    ];             
problem.dmax =  [   inf;    inf;    0.44;   0.6;    inf;    inf];
problem.dmin = -[   inf;    inf;    0.44;   0.6;    inf;    inf];

problem.Q  = eye(6);           
problem.Q(1,1) = 100;
problem.Q(2,2) = 100;
problem.Q(3,3) = 10;
problem.Q(4,4) = 10;
problem.Q(5,5) = 400;
problem.Q(6,6) = 200;
problem.R  = 0.001*eye(2);
problem.S = zeros(6,2);
problem.C = eye(6);

[~, problem.P] = dlqr(problem.A, problem.B, problem.Q, ...
    problem.R, problem.S);     % optional (default, zero matrix)

problem.x0 = [0.5,0.5,0,0,0,0]' ;    

problem.ni   = 20;                
problem.ur   = {[0;0]};
problem.yr = {[0;0;0;0;0;0]};
problem.xNr = [0;0;0;0;0;0];

end
