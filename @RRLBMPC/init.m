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
    obj = obj.getHmatrix;
   
end

function obj = getHmatrix(obj)
    % get the H1 and H2 matrices for the ALADIN method
    Q_bar = obj.Q;
    % The Hessian matrix of the function b_tilde
    Hessian_b = diag(obj.Wx) / obj.delta_relax;
    
    Q_bar = Q_bar + obj.rho * Hessian_b;

    obj.ALADIN_H1 = kron(eye(obj.N), Q_bar);
    obj.ALADIN_H2 = kron(eye(obj.N), obj.R);

    

end