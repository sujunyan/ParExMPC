
 function [z1_next, z2_next, lam_next, u0] = ADMM_solve(obj, z1, z2, lam, maxiter)
    % Run multiple iterations of ADMM.

    for iter = 1:maxiter
        [z1_next, z2_next, lam_next, u0] = ADMM_one_iteration(obj, z1, z2, lam);
        % fprintf("ADMM iteration %d completed. z1 error %.2f z2 error %.2f lam_diff %.2f \n", iter, norm(z1-z1_next), norm(z2 - z2_next), norm(lam - lam_next));
        z1 = z1_next;
        z2 = z2_next;
        lam = lam_next;
    end
end

function res = f1_ADMM(obj, z1, z2, lam, sigma)
    % The objective of the first subproblem in ADMM
    a1 = f1(obj, z1);
    a2 = lam' * z1;
    a3 = ( sigma / 2 ) * norm( z1 - obj.compactA * z2 - obj.compactb )^2;
    res = a1 + a2 + a3;
end

function [z1_next, z2_next, lam_next, u0] = ADMM_one_iteration(obj, z1, z2, lam)
    sigma = obj.ADMM_sigma;

    % Solve the first subproblem of ADMM
            
    options = optimset('Display', 'off');
    [z1_next, fval1] = fminunc(@(z1_var) f1_ADMM(obj, z1_var, z2, lam, sigma), z1, options);
    % fval1

    C = kron(speye(obj.N), obj.Cu);
    d = kron(ones(obj.N,1), obj.du);
    z2_next = fmincon(@(z2_var) f2_ADMM(obj, z1_next, z2_var, lam, sigma), z2, C, d, [], [], [], [], [], options);


    lam_next = lam + sigma * ( z1_next - obj.compactA * z2_next - obj.compactb );

    u0 = getZ2k(obj, z2_next, 1);

end

function res = f2_ADMM(obj, z1, z2, lam, sigma)
    % The objective of the second subproblem in ADMM
    a1 = f2_no_cons(obj, z2) - (obj.compactA * z2)' * lam;
    a2 = ( sigma / 2 ) * norm( z1 - obj.compactA * z2 - obj.compactb)^2;
    res = a1 + a2;
end
