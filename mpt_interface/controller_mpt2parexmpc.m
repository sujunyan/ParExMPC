function [ mpc_ctrl ] = controller_mpt2parexmpc( mpt_ctrl, cons_mul, mptSolver )
%% CONTROLLER_MPT2PAREXMPC
%
% Function CONTROLLER_MPT2PAREXMPC translates MPC controller constructed in
% MPT toolbox into ParExMPC controller. 
%
% [ mpc_ctrl ] = controller_mpt2parexmpc( mpt_ctrl, cons_mul, mptSolver )
%
% where:
%
% mpt_ctrl - is input MPC controller constructed in MPT toolbox
% N        - is input scalar defining the prediction horizon 
% mpc_ctrl - is output object of ParExMPC defining the MPC controller
%
% 2023-12-03

% Default inputs
if( nargin < 3 )
    mptSolver = 'plcp';
end
if( nargin < 2 )
    cons_mul = '1';
end

mpc_ctrl = peMPC(mpt_ctrl.model.A, mpt_ctrl.model.B, ...
    mpt_ctrl.model.x.penalty.H , mpt_ctrl.model.u.penalty.H, ... 
    mpt_ctrl.model.x.terminalPenalty.H, ...
    'umin', mpt_ctrl.model.u.min, 'umax', mpt_ctrl.model.u.max, ...
    'dmin', mpt_ctrl.model.x.min * cons_mul, 'dmax', mpt_ctrl.model.x.max * cons_mul, ... % TODO: double-check state constraints 
    'N', mpt_ctrl.N, ...
    'cons_mul', cons_mul, 'mptSolver', mptSolver); % Optional ipnuts

end