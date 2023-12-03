function [ mpc_ctrl ] = mpt2parexmpc( obj, varargin )
%% MPT2PAREXMPC
%
% Function MPT2PAREXMPC translates LTI system model formulated in MPT
% toolbox into ParExMPC controller.
%
% [ mpc_ctrl ] = mpt2parexmpc( mpt_ctrl )
% [ mpc_ctrl ] = mpt2parexmpc( model, 'N', N )
% [ mpc_ctrl ] = mpt2parexmpc( mpt_ctrl/model, 'N', N, 'cons_mul', cons_mul, 'mptSolver', mptSolver )
%
% where:
%
% mpt_ctrl - is input MPC controller constructed in MPT toolbox
% model     - is input LTI system defined in MPT toolbox
% N         - is input scalar defining the prediction horizon
% cons_mul  - is optional input scalar scaling (state) constraints (default: 1)
% mptSolver - is optional input string of multiparametric solver (default: 'plcp')
% mpc_ctrl  - is output object of ParExMPC defining the MPC controller
%
% 2023-12-03

%% Input parser
% TODO: More arguments
p = inputParser;
% addOptional(p,'xr',zeros(obj.nx,1));
% addOptional(p,'xNr',zeros(obj.nx,1));
% addOptional(p,'ur',zeros(obj.nu,1));
% addOptional(p,'C',eye(obj.nx));
% addOptional(p,'D',zeros(obj.nx,obj.nu));
% addOptional(p,'dmin',-inf*ones(obj.nx,1));
% addOptional(p,'dmax',inf*ones(obj.nx,1));
% addOptional(p,'umin',-inf*ones(obj.nu,1));
% addOptional(p,'umax',inf*ones(obj.nu,1));
addOptional(p,'N', -Inf);
addOptional(p,'cons_mul',1);
addOptional(p,'mptSolver','plcp');

parse(p,varargin{:});
% obj.xr = p.Results.xr;
% obj.xNr = p.Results.xNr;
% obj.ur = p.Results.ur;
% obj.C = p.Results.C;
% obj.D = p.Results.D;
% umin = p.Results.umin;
% umax = p.Results.umax;
% dmin = p.Results.dmin;
% dmax = p.Results.dmax;
N = p.Results.N;
cons_mul = p.Results.cons_mul;
% dmin = dmin * cons_mul;
% dmax = dmax * cons_mul;
mptSolver = p.Results.mptSolver;

MyClass = class( obj );

if( isequal( MyClass(end-8:end),'ontroller' ) == 1 )
    [ mpc_ctrl ] = controller_mpt2parexmpc( obj, cons_mul, mptSolver );
elseif( isequal( MyClass(end-8:end),'LTISystem' ) == 1 )
    [ mpc_ctrl ] = model_mpt2parexmpc( obj, N, cons_mul, mptSolver  );
else
    error(sprintf('MPT2PAREXMPC: Unexpected class of input object: "%s"!',MyClass))
end

end