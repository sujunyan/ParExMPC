function obj = solveDecouple(obj,x0)
% solve the decoupled problems.
% Output:
%   xi: 1xN cell that store the solution of the decompled problems
N = obj.N;

xiN = solveDecoupleN(obj,x0);
xi{N+1} = xiN;
Hs = 4*obj.Sigmak;
%As =  [obj.C*[obj.A,obj.B]; blkdiag(obj.C,eye(obj.nu))];
As = [obj.C,obj.D; zeros(obj.nu,obj.nx), eye(obj.nu)];
lbs = [];
ubs = [];
%lbAs = [obj.xmin;obj.xmin;obj.umin];
%ubAs = [obj.xmax;obj.xmax;obj.umax];
lbAs = [obj.dmin;obj.umin];
ubAs = [obj.dmax;obj.umax];
A = [As;-As];
b = [ubAs;-lbAs];
quadprog_opt = optimoptions('quadprog','Display','off');
for kk = 1:obj.N
    gs = [-2*obj.Q*obj.xr;-2*obj.R*obj.ur] -0.5*Hs*obj.z{kk} ...
        - (obj.Gk'*obj.lam{kk} - obj.Hk'*obj.lam{kk+1}); % The sign is inconsistent with the paper
    xi{kk} = quadprog(Hs,gs,A,b,[],[],[],[],[],quadprog_opt);
end
obj.xi = xi;
end % end of solveDecouple

function xi = solveDecoupleN(obj,x0)
% solve the Last decoupled problem with quadprog
Hs = 4*obj.SigmaN;
%gs = -0.5*Hs*(obj.xNr + obj.z{obj.N+1}) - (obj.GN'*obj.lam{obj.N});% The sign is inconsistent with the paper
gs = -2*obj.P*(obj.xNr + obj.z{obj.N+1}) - obj.lam{obj.N+1};% The sign is inconsistent with the paper
As = obj.C;
lbs = [];
ubs = [];
lbAs = obj.dmin;
ubAs = obj.dmax;
A = [As;-As];
b = [ubAs;-lbAs];
quadprog_opt = optimoptions('quadprog','Display','off');
xi = quadprog(Hs,gs,A,b,[],[],[],[],[],quadprog_opt);
end