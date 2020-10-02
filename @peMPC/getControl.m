function [obj,u,x] = getControl(obj,x0,z,lam,max_iter)
    % perform the aladin algorithm to the MPC optimization problem without mex functions
    %
    % Input:
    %   x0: The initial state 
    %   z:  initial guess of primal variable
    %   lam: initial guess of dual variable
    % Output:
    %   u0: The control input computed by peMPC controller

    % TODO: should rescale first but skip it for now.
    if nargin >= 4 % add initial guess 
        % warm start
        lam = reshape(lam,obj.nx,[]);
        for kk = 1:obj.N
            obj.lam{kk} = lam(:,kk);
        end
        y1_N = reshape(z(1:(obj.nx+obj.nu)*obj.N),obj.nx+obj.nu,[]);
        for kk = 1:obj.N
            obj.z{kk} = y1_N(:,kk);
        end
        obj.z{obj.N+1} = z((obj.nx+obj.nu)*obj.N+1:end);
    end
    if nargin < 5
        max_iter = obj.maxiter;
    end


    obj.KKT_con_res = [-obj.A*x0; zeros((obj.nx)*(obj.N-1),1)];
    for m = 1:max_iter
        obj = obj.solveDecouple(x0);
        obj = obj.solveCoupleRiccati(x0);
        if m ~= max_iter
            xi_v = vertcat(obj.xi{:});
            z_v = vertcat(obj.z{:});
            gap = norm(xi_v-z_v,inf);
            %fprintf("The gap is %e\n",gap);
            if gap < 1e-3
                break
            end
        end
    end

    %% construct the output ---------------------------------------
    u0 = obj.xi{1};
    u0 = u0(obj.nx+1:end);
    if nargout == 2
        u = u0;
    else
        xi_v = vertcat(xi{:});
        ny = obj.nu+obj.nx;
        uxz = xi_v(obj.nu+1:obj.nu+ny*(obj.N-1));
        uxz = reshape(uxz,ny,[]);
        xN = xi_v(end-obj.nx+1:end);
        x = [x0,uxz(1:obj.nx,:),xN];
        u = [u0,uxz(obj.nx+1:obj.nx+obj.nu,:)];
    end
    % shift the variable such for next stage --------
    obj.z = {obj.z{2:end-1},[obj.z{end};zeros(obj.nu,1)],zeros(obj.nx,1)};
    obj.lam = {obj.lam{2:end},zeros(obj.nx,1)};
end