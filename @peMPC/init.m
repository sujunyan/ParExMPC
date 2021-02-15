function obj = init(obj)
  % initilize  
  obj.nx = size(obj.B,1);
  obj.nu = size(obj.B,2);
  obj = obj.getGH;
  obj = obj.getGamma;
  obj = obj.getMaxIter;
  obj.isFirst = true;
  %obj.xr = obj.xr'*obj.C ;
  obj.zr = [kron(ones(obj.N,1),[obj.xr;obj.ur]) ;obj.xNr];
  tol = 1e-8; % TODO
  %obj.Q = obj.C'*obj.Q*obj.C; % for low rank C, we first convert it to nonnegative square matrix
  obj.R = obj.R + tol*eye(size(obj.R)); 
  obj.Q = obj.Q + tol*eye(size(obj.Q)); 
  obj.P = obj.P + tol*eye(size(obj.P));
  obj.Sigma0 = obj.R;
  obj.Sigmak = blkdiag(obj.Q,obj.R);
  obj.Sigmak = obj.Sigmak + tol*eye(size(obj.Sigmak));
  obj.NSigmak = kron(speye(obj.N-1),obj.Sigmak);
  obj.SigmaN = obj.P;
  % TODO: add robustness to avoid unbound decoupled problems
  % init the variables
      
  obj.z{obj.N+1} = zeros(obj.nx,1);
  obj.lam{obj.N+1} = zeros(obj.nx,1);
  for kk = 1:obj.N
      obj.lam{kk} = zeros(obj.nx,1);
      obj.z{kk} = zeros(obj.nu+obj.nx,1);
  end
  %obj = obj.preCQP;
  obj = obj.preRiccati;
  obj = obj.getLargeQP;
end