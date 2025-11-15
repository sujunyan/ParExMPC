
 function [z1_next, z2_next, lam_next, u0] = ALADIN_solve(obj, z1, z2, lam, maxiter)
    % Run multiple iterations of ADMM.

    for iter = 1:maxiter
        [z1_next, z2_next, lam_next, u0] = ALADIN_one_iteration(obj, z1, z2, lam);
        z1 = z1_next;
        z2 = z2_next;
        lam = lam_next;
    end
end


function [z1_next, z2_next, lam_next, u0] = ALADIN_one_iteration(obj, z1, z2, lam)
    H1 = obj.ALADIN_H1; H2 = obj.ALADIN_H2;

    options = optimset('Display', 'off');
    [xi1, fval1] = fminunc(@(xi1_var) f1_ALADIN(obj, xi1_var, z1, lam), z1, options);

    g1 = obj.ALADIN_H1 * (z1 - xi1) - lam;


    C = kron(speye(obj.N), obj.Cu);
    d = kron(ones(obj.N,1), obj.du);
    [xi2, fval2] = fmincon(@(xi2_var) f2_ALADIN(obj, xi2_var, z2, lam), z2, C, d, [], [], [], [], [], options);

    g2 = obj.ALADIN_H2 * (z2 - xi2) + obj.compactA' * lam;

    % solve for the coupled QP
    H = blkdiag(obj.ALADIN_H1, obj.ALADIN_H2);

    g11 = g1 - xi1' * H1;
    g22 = g2 - xi2' * H2;
    g = vertcat(g11, g22);

    Aeq = [eye(length(z1)), -obj.compactA];
    beq = obj.compactb;
    z0 = vertcat(z1, z2);

    [z, fval3, exitflag, output, lambda] = quadprog(H, g, [], [], Aeq, beq, [], [], z0, options);

    z1_next = z(1:length(z1));
    z2_next = z(length(z1)+1:end);
    lam_next = lambda.eqlin;
    u0 = getZ2k(obj, z2_next, 1);
end

function res = f1_ALADIN(obj, xi1, z1, lam)
    % The objective of the first subproblem in ADMM
    H1 = obj.ALADIN_H1;
    a1 = f1(obj, xi);
    a2 = lam' * z1;
    xi_diff = xi1 - z1;
    a3 = 1/2 * xi_diff' * H1 * xi_diff;
    res = a1 + a2 + a3;
end

function res = f2_ALADIN(obj, xi2, z2, lam)
    a1 = f2_no_cons(obj, xi2);
    a2 = - (obj.compactA * xi2)' * lam;
    a3 = 1/2 * (xi2 - z2)' * obj.ALADIN_H2 * (xi2 - z2);
    res = a1 + a2 + a3;
end

function obj = solveCoupleUpdate(obj,x0)
    % Solve the coupled QP problem of the following form
    % and update the primal and dual variable
    % min 1/2 x'Qx + q'x; s.t. Ax=b | lam and the corresponding KKT matrix is given by
    % [Q,A'; * [x;   = [-q;
    %  A,0]     lam]    b]
    ny = obj.nu+obj.nx;
    xi_v = vertcat(obj.xi{:}); % the vectorized xi
    z_v = vertcat(obj.z{:});
    [del_lam,primal] = pempc_solve_couple_mex(xi_v,z_v,x0);
    %tmp = 2*xi_v - y_v; %+ obj.yr; % TODO
    %obj.KKT_res = [obj.KKT_obj*tmp;obj.KKT_con_res];
    %sol_cQP = obj.dKKT\obj.KKT_res;
    %primal = sol_cQP(1:ny*obj.N);
    %del_lam  = sol_cQP(ny*obj.N+1:end);
    %obj.z{1} = primal(1:obj.nu);
    obj.z{obj.N+1} = primal(end-obj.nx+1:end);
    for kk = 1:obj.N
        istart = (kk-1)*(ny)+1;
        obj.z{kk} = primal(istart:istart+ny-1);
    end
    nxz = obj.nx;
    for kk = 1:obj.N+1
        istart = (kk-1)*(nxz)+1;
        obj.del_lam{kk} = del_lam(istart:istart+nxz-1);
        obj.lam{kk} = obj.lam{kk} + obj.del_lam{kk};
    end
end