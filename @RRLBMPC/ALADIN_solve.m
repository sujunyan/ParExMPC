
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

   

end