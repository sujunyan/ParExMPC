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

    obj.lam = zeros((obj.N+1)*obj.nx,1);
   
end