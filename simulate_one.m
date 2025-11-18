
function res = simulate_one(mpc0, nsim, sim_tol, method, max_iter0, x0)
    % simulate a single run of the MPC problem with given parameters


    nx = mpc0.nx; nu = mpc0.nu; N = mpc0.N;
    time = 0; % time used by QP solver
    x_vec = [x0];
    u0_vec = [];
    J_vec = [0];
    for i = 1:nsim
        fprintf(" Simulation: %d/%d \n", i, nsim);
        x0 =  x_vec(:,end);
        mpc0 = mpc0.updateX0(x0);
        max_iter = max_iter0;
        if i == 1
            z1 = kron(ones(N,1), x0);
            z2 = zeros(nu * N, 1);
            lam = zeros(nx * N, 1);
            max_iter = 100;
        end
        tic;
        if method == "ADMM"
            [z1, z2, lam, u0] = mpc0.ADMM_solve(z1, z2, lam, max_iter);
        elseif method == "ALADIN"
            [z1, z2, lam, u0] = mpc0.ALADIN_solve(z1, z2, lam, max_iter);
        elseif method == "fmincon"
            [z1, z2, u0] = mpc0.fmincon_solve(z1, z2);
        end
        elapsed = toc;
        time = time + elapsed;
        xn = mpc0.A*x0 + mpc0.B*u0;
        xn_sim = z1(1:nx);
        xn ;
        z2 ;
        
        
        x_vec = [x_vec, xn];
        u0_vec = [u0_vec, u0];
        Jn = J_vec(end) + x0'*mpc0.Q*x0 + u0'*mpc0.R*u0;
        J_vec = [J_vec, Jn];
        if norm(xn - x0, Inf) < sim_tol
            fprintf(" Converged at step %d \n", i);
            break;
        end
    end

    res = struct("x_vec", x_vec, "J_vec", J_vec, "time", time, "u0_vec", u0_vec, "isim_stop", i);

end