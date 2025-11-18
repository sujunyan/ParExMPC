
function [z1_next, z2_next, u0] = fmincon_solve(obj, z1, z2)
    % Directly solve the MPC problem with the fmincon solver

    nx = obj.nx; mu = obj.mu; N = obj.N;
    options = optimset('Display', 'off');

    z0 = vertcat(z1, z2);

    Aeq = [eye(obj.nx * obj.N), -obj.compactA];
    beq = obj.compactb;
    A = [zeros(nx*N, mu*N) , kron(speye(obj.N), obj.Cu)];
    b = kron(ones(obj.N,1), obj.du);

    z = fmincon(@(z_var) fmincon_func(obj, z_var), z0, A, b, Aeq, beq, [], [], [], options);

    z1_next = z(1:obj.nx * obj.N);
    z2_next = z(obj.nx * obj.N + 1:end);
    u0 = getZ2k(obj, z2_next, 1);
    
end

function res = fmincon_func(obj, z)
    z1_length = obj.nx * obj.N;
    z1 = z(1:z1_length);

    z2 = z(z1_length+1:end);

    a1 = f1(obj, z1);
    a2 = f2_no_cons(obj, z2);

    res = a1 + a2;
end