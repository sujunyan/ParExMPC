function obj = init(obj)
    % initilize
    obj.nx = size(obj.B,1);
    obj.nu = size(obj.B,2);

    obj = obj.getWx;
    obj.isFirst = true;

    tol = 1e-8; % TODO
    %obj.Q = obj.C'*obj.Q*obj.C; % for low rank C, we first convert it to nonnegative square matrix
    obj.R = obj.R + tol*eye(size(obj.R));
    obj.Q = obj.Q + tol*eye(size(obj.Q));
    obj.P = obj.P + tol*eye(size(obj.P));

    obj = obj.getCompactForm;
    % obj.lam = zeros((obj.N+1)*obj.nx,1);
    % obj = obj.getHmatrix;
    obj = getHmatrix(obj);
   
end

function obj = getHmatrix(obj)
    % get the H1 and H2 matrices for the ALADIN method
    Q_bar = obj.Q;
    % The Hessian matrix of the function b_tilde
    % Hessian_b = diag(obj.wx) / obj.delta_relax;
    Hessian_b = zeros(obj.nx, obj.nx);
    for r = 1:obj.mx
        Cr = obj.Cx(r,:);
        dr = obj.dx(r);
        wr = obj.wx(r);
        H_tmp = Cr * Cr' * obj.relax_barrier_d2(obj.delta_relax, dr);
        Hessian_b = Hessian_b + H_tmp * wr * obj.rho;
    end
    
    Q_bar = Q_bar + obj.rho * Hessian_b;

    obj.ALADIN_H1 = kron(eye(obj.N), Q_bar);
    obj.ALADIN_H2 = kron(eye(obj.N), obj.R);

    

end