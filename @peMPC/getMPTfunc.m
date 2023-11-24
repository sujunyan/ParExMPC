function obj = getMPTfunc(obj)
% get MPT functions for decouple problems

% The bound of the parameter theata
fprintf("Calling MPT functions...\n");

mptopt('verbose',0);
global MPTOPTIONS;
max_infbound = 5e0*norm(blkdiag(obj.Sigmak,obj.SigmaN)) ...
    *max([norm(obj.H0),norm(obj.Hk),norm(obj.Gk),norm(obj.GN)])...
    * norm( [obj.zmax(isfinite(obj.zmax)) ] ,inf) ;
obj.MPT_inf_bound = max_infbound;
%fprintf("In getMPTfunc, get max_infbound %e\n",max_infbound);
mptopt('infbound',max_infbound);
mptopt('pqpsolver', obj.mptSolver); % list of solver: plcp, mpqp, enumplcp, enumpqp, rlenumpqp
class_folder = fileparts(mfilename('fullpath'));

tiebreak = 'obj';
tstart = tic;
% comment out the x0 since we dualize the x0
% for lam0 ------------------------------------------------------
%xi0 = sdpvar(obj.nu,1);
%gs = sdpvar(obj.nu,1);
%x0 = sdpvar(obj.nx,1);
%Hs = 4*obj.Sigma0;

%pF = [zeros(obj.nu,obj.nx), eye(obj.nu)];
%% A*u <= b + pB*th
%b = [obj.xmax;-obj.xmin];
%pB = [-obj.A,zeros(obj.nx,obj.nu);obj.A,zeros(obj.nx,obj.nu)];
%mpt0A = [obj.B;-obj.B];
%% add upper/lower bound of parameters
%%mpt0A = [mpt0A;zeros(2*(obj.nu+obj.nx),obj.nu)];
%%pB = [pB; -eye(obj.nu+obj.nx); eye(obj.nu+obj.nx)];
%%b = [b; max_infbound*ones(obj.nu+obj.nx,1); -max_infbound*ones(obj.nu+obj.nx,1)];
%pB = pB(isfinite(b),:);
%mpt0A = mpt0A(isfinite(b),:);
%b = b(isfinite(b),:);
%
%problem0 = Opt('H',Hs,'pF',pF,'lb',obj.umin,'ub',obj.umax,...
%            'A',mpt0A,'b',b,'pB',pB);
%res0 = problem0.solve;
%res0.xopt.toMatlab([func_folder filesep 'MPT_func0.m'],'primal',tiebreak);

% for lamN --------------------------------------

mptopt('pqpsolver', obj.mptSolver); % list of solver: plcp, mpqp, enumplcp, enumpqp, rlenumpqp
lamN = sdpvar(obj.nx,1);
yN = sdpvar(obj.nx,1);
xiN = sdpvar(obj.nx,1);
Hs = 4*obj.SigmaN;
gs = sdpvar(obj.nx,1);
%gs = -0.5*Hs*(2*obj.xNr + yN) - (obj.GN'*lamN);% The sign is inconsistent with the paper
lbs = obj.dmin;
ubs = obj.dmax;
gs_max = MPTOPTIONS.infbound*ones(size(gs));
gs_min = -gs_max;
objective = 1/2*xiN'*Hs*xiN + gs'*xiN;
constraints = [gs_min<=gs<=gs_max, lbs <= obj.C*xiN <= ubs ];
problemN = Opt(constraints,objective,gs,xiN);
%problemN = Opt('H',Hs,'pF',eye(obj.nx),'lb',obj.xmin,'ub',obj.xmax);
fprintf("Calling MPT functions for N-th problem -----------------------------------------------------\n");
resN = problemN.solve;
%resN.xopt.toMatlab([func_folder filesep 'MPT_func_N.m'],'primal',tiebreak);

% for lamk -------------------------------------------------------------
mptopt('pqpsolver', obj.mptSolver); % list of solver: plcp, mpqp, enumplcp, enumpqp, rlenumpqp

xik = sdpvar(obj.nx+obj.nu,1);
Hs = 4*obj.Sigmak;
gs = sdpvar(obj.nx+obj.nu,1);
gs_max = MPTOPTIONS.infbound*ones(size(gs));
gs_min = -gs_max;
%gs = -0.5*Hs*(2*[obj.xr;obj.ur;obj.zr] + yk) ...
%- (obj.Gk'*lamk_1 - obj.Hk'*lamk); % The sign is inconsistent with the paper

As = [obj.C,obj.D; zeros(obj.nu,obj.nx), eye(obj.nu)];
lbAs = [obj.dmin;obj.umin];
ubAs = [obj.dmax;obj.umax];

objective = 1/2*xik'*Hs*xik + gs'*xik;
constraints = [gs_min<= gs<= gs_max, lbAs<= As*xik <=ubAs];
%constraints = [lbAs<= As*xik <=ubAs];
problemk = Opt(constraints,objective,gs,xik);
fprintf("Calling MPT functions for k-th problem -----------------------------------------------------\n");
resk = problemk.solve;

fprintf("Multiparametric problem solved in %f seconds.\n",toc(tstart));
fprintf("Export to C-code:\n");

fprintf("Export of file 'mpQP_data.h'.\n");
filewrite_t = tic;
fileID = fopen([class_folder filesep 'c_code' filesep 'mpQP_data.h'],'w+');
fprintf(fileID,"/*predefined mpt matrices for explict evaluation*/\n");
fprintf(fileID,"#ifndef __PRE_MPT_H__\n#define __PRE_MPT_H__\n");
%obj.PolyUnionToC(fileID, {res0.xopt},'z');
obj.PolyUnionToC(fileID, {resN.xopt},'N');
obj.PolyUnionToC(fileID, {resk.xopt},'k');
fprintf(fileID,"#endif\n");
fclose(fileID);
fprintf("Export to 'pre_mpt.h' required %f seconds.\n",toc(filewrite_t));
fprintf("Calling MPT functions required %f seconds.\n",toc(tstart));
end