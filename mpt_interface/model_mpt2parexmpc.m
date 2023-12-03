function [ mpc_ctrl ] = model_mpt2parexmpc( model, N, cons_mul, mptSolver )
%% MODEL_MPT2PAREXMPC
%
% Function MODEL_MPT2PAREXMPC translates LTI system model formulated in MPT
% toolbox into ParExMPC controller.
%
% [ mpc_ctrl ] = model_mpt2parexmpc( model, N, cons_mul, mptSolver )
%
% where:
%
% model    - is input LTI system defined in MPT toolbox
% N        - is input scalar defining the prediction horizon 
% mpc_ctrl - is output object of ParExMPC defining the MPC controller
%
% 2023-12-03

% Default inputs
if( nargin < 4 )
    mptSolver = 'plcp';
end
if( nargin < 3 )
    cons_mul = '1';
end

mpc_ctrl = peMPC(model.A, model.B, ...
    model.x.penalty.H , model.u.penalty.H, model.x.terminalPenalty.H, ....
    'umin', model.u.min, 'umax', model.u.max, ...
    'dmin', model.x.min * cons_mul, 'dmax', model.x.max * cons_mul, ... % TODO: double-check state constraints 
    'N', N, ...
    'cons_mul', cons_mul, 'mptSolver', mptSolver); % Optional ipnuts

end