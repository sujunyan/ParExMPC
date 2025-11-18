
function [z1_next, z2_next, u0] = quadprog_state_solve(obj, z1, z2)
    % Directly solve the MPC problem with the fmincon solver
    % In this function varaint, we also consider both state and input constraints as hard constraints
    % This is to compare the impact from relaxed state constraint

    nx = obj.nx; mu = obj.mu; N = obj.N;
    options = optimset('Display', 'off');

    z0 = vertcat(z1, z2);

    Aeq = [eye(obj.nx * obj.N), -obj.compactA];
    beq = obj.compactb;
    % A = [zeros(nx*N, mu*N) , kron(speye(obj.N), obj.Cu)];
    A = blkdiag(  kron(speye(N), Cx),  kron(speye(N), Cu))
    b1 = kron(ones(obj.N,1), obj.dx);
    b2 = kron(ones(obj.N,1), obj.du);
    b = vertcat(b1, b2);

    Q = blkdiag(obj.compactQ, obj.compactR);
    q = zeros(size(z0));
    

    z = quadprog(Q, q, A, b, Aeq, beq, [], [], z0, options)

    z1_next = z(1:obj.nx * obj.N);
    z2_next = z(obj.nx * obj.N + 1:end);
    u0 = getZ2k(obj, z2_next, 1);
    
end

