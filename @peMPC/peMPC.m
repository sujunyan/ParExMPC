classdef peMPC
    % parallel explicit MPC controller
    % TODO: refer to the paper: https://arxiv.org/abs/1903.06790
    % Author:         Junyan Su
    % The pempc object to solve the problem of the form
    %   \min \sum (x[k]-x_r)' Q (x[k]-x_r) for k = 1...N-1
    %           + (u[k]-u_r)' R (u[k]-u_r)
    %           + (x[N]-x_Nr)' P (x[N]-x_Nr)
    %   s.t. x[k+1] = Ax[k] + Bu[k]
    %        dmin  <= Cx[k] + Du[k] <= dmax
    %        umin  <= u[k] <= umax
    % Optional input:
    %       xr: the reference state
    %       xNr: the reference terminal state
    %       ur: the reference control input
    %       dmin,dmax: the state contraint
    %       umin,umax: the control input contraint
    %       N: the time horizon
    %       cons_mul: the contraint multipiler, often set less than 1 to avoid contraint violation
    %       par_flag: the boolean flag to enable/disable the parallism implementaion
    %       par_threshold: if the time horizon is larger than this threshold, then enable the parallel computing, valid only when par_flag is true.
    % Call:
    %       mpc0 = peMPC(A,B,Q,R,P,{,optional inputs})

    properties
        % The System dynamics
        % x[k+1] = A x[k] + B u[k]
        % dmin <= C*x + D*u <= dmax
        A
        B
        C
        D
        % TODO: for now, we assume the system is fully observable

        % The converted system dynamics---in case we have a different system format
        % G[k+1]z[k+1] = H[k]z[k] + h[k]
        % where z[0] = u[0]; z[k] = [x[k]; u[k]]; z[N] = x[N]
        % GN
        Gk
        NGkt        % kron(eye(N-1),Gk')
        GN
        Hk
        NHkt        % kron(eye(N-1),Hk')
        H0

        % Objective function
        Q           %  (x[k]-x_r[k])' Q (x[k]-x_r[k]) for k = 1...N-1
        R           %  (u[k]-u_r[k])' R (u[k]-u_r[k])
        P           %  (x[N]-x_r[N])' P (x[N]-x_r[N])
        Sigma0      % blkdiag(R,S)
        Sigmak      % blkdiag(Q,R,S)
        NSigmak     % kron(eye(N-1),Sigmak)
        SigmaN      % P

        % Reference values
        % TODO: for now, we treat it as constant
        xr          % reference state trajectory
        ur          % reference input trajectory
        xNr         % reference terminal state
        zr          % The reference stack variable, have size (nx+nu)*N

        % State and input constraints
        dmin        % dmin <= C*x + Du <= dmax
        dmax
        umin        % umin <= u[k] <= umax
        umax
        z           % The stacked variable [x;u]
        zmin
        zmax

        % the dimentions
        nx
        nu

        % controller parameters
        maxiter     % the maximum iteration
        N           % the time horizon
        gamma       % the rescale factor
        isFirst     % the flag to indicate that if this object has ben used. For the first time, we run a large number of itertions in aladin to initilize the peMPC.

        % the stored variable
        lam         % the Lagrangian multiplier
        del_lam     % TODO: for testing
        xi          % the alternating direction

        % precomputed coupled QP KKT matrix
        dKKT        % The decomposed KKT matrix for faster speed
        KKT
        KKT_obj
        KKT_con_res
        KKT_con
        KKT_res
        % The variables used in Riccati based method
        Ric_P
        Ric_Lam
        Ric_Lam_inv
        Ric_L

        % properties for solving the qp in condensing form
        % for comparison and inilization
        % min 1/2 x'Hx + g'x
        % s.t. QP_l <= Ax <= QP_u
        QP_H
        QP_g
        QP_A
        QP_l
        QP_u

        % additional variables -----------------------------------
        MPT_P       % the polyhedra Union solved by MPT, used in the parallel step.
        MPT_inf_bound % the infinite bound of MPT
        % The MPT solver to be used
        % list of solver: plcp, mpqp, enumplcp, enumpqp, rlenumpqp
        % default is plcp, try enumpqp for better speed (might fail for some cases)
        mptSolver
        cons_mul    % constraint multipiler to shrink the constraints and avoid constraint violation
        use_parallel % a flag to choose if we want to use parallelism or not.
        parallel_threshold % if the time horizon larger than this threshold, then enable the parallel computing.

    end % End of the properties

    methods (Access = public)

        function out = isPosDef(obj,A)
            % determine if matrix A is positive definite.
            tol = 1e-8;
            d = eig(A);
            out = all(d > tol);
        end

        function obj = peMPC(A,B,Q,R,P,varargin)
            % The constructor of the peMPC
            % get from a problem object
            % input:
            %   A,B: The system dynamics x[k+1] = Ax[k]+Bu[k]
            %   Q,R,P: The objective function
            %       (x[k]-x_r)' Q (x[k]-x_r) for k = 1...N-1
            %       (u[k]-u_r)' R (u[k]-u_r)
            %       (x[N]-x_Nr)' P (x[N]-x_Nr)
            % Optional input:
            %   xr: the reference state
            %   xNr: the reference terminal state
            %   ur: the reference control input
            %   xmin,xmax: the state contraint
            %   umin,umax: the control input contraint
            %   N: the time horizon
            %   cons_mul: the contraint multipiler, often set less than 1 to avoid contraint violation
            %   par_flag: the boolean flag to enable/disable the parallism
            %   implementaion 
            %   
            % Call:
            %   mpc0 = peMPC(A,B,Q,R,P,{,optional inputs})
            if (nargin == 0)
                fprintf("Empty peMPC object created\n");
                return;
            end
            obj.A = A; obj.B = B; obj.Q = Q; obj.R = R; obj.P = P;
            if ~obj.isPosDef(obj.Q)
                warning("The matrix Q is not strictly positive definite, the result might be wrong");
            end
            if ~obj.isPosDef(obj.P)
                warning("The matrix P is not strictly positive definite, the result might be wrong");
            end
            if ~obj.isPosDef(obj.R)
                warning("The matrix R is not strictly positive definite, the result might be wrong");
            end

            obj.nx = size(obj.A,1);
            obj.nu = size(obj.B,2);
            % TODO: More arguments
            p = inputParser;
            addOptional(p,'xr',zeros(obj.nx,1));
            addOptional(p,'xNr',zeros(obj.nx,1));
            addOptional(p,'ur',zeros(obj.nu,1));
            addOptional(p,'C',eye(obj.nx));
            addOptional(p,'D',zeros(obj.nx,obj.nu));
            addOptional(p,'dmin',-inf*ones(obj.nx,1));
            addOptional(p,'dmax',inf*ones(obj.nx,1));
            addOptional(p,'umin',-inf*ones(obj.nu,1));
            addOptional(p,'umax',inf*ones(obj.nu,1));
            addOptional(p,'N',10);
            addOptional(p,'cons_mul',1);
            addOptional(p,'mptSolver','plcp');
            addOptional(p,'par_flag', true);
            addOptional(p,'par_threshold', 20);

            parse(p,varargin{:});
            obj.xr = p.Results.xr;
            obj.xNr = p.Results.xNr;
            obj.ur = p.Results.ur;
            obj.C = p.Results.C;
            obj.D = p.Results.D;
            obj.umin = p.Results.umin;
            obj.umax = p.Results.umax;
            obj.dmin = p.Results.dmin;
            obj.dmax = p.Results.dmax;
            obj.N = p.Results.N;
            obj.cons_mul = p.Results.cons_mul;
            obj.dmin = obj.dmin * obj.cons_mul;
            obj.dmax = obj.dmax * obj.cons_mul;
            obj.mptSolver = p.Results.mptSolver;
            obj.use_parallel = p.Results.par_flag;
            obj.parallel_threshold = p.Results.par_threshold;
        end

        function obj = getGH(obj)
            % convert the system dynamics to G and H
            nx = obj.nx;
            nu = obj.nu;
            A = obj.A; B = obj.B;
            obj.H0 = [B];
            obj.Hk = [A,B];
            obj.Gk = [eye(nx), zeros(nx,nu)];
            obj.GN = eye(nx); % TODO: inconsistent with the paper
            % the constraint for the augmented variable
            obj.NGkt = kron(speye(obj.N-1),obj.Gk');
            obj.NHkt = kron(speye(obj.N-1),obj.Hk');
            obj.zmin = [obj.dmin; obj.umin];
            obj.zmax = [obj.dmax; obj.umax];
        end

        function obj = getGamma(obj)
            % TODO
        end

        function obj = getMaxIter(obj, varargin)
            % Function "getMaxIter" evaluates minimum necessary iterations
            % for the stability guarantees (m_bar)
            %
            % obj = getMaxIter(obj, kappa, gamma, sigma, eta, tau)
            % 
            % maxiter = 2*log( 2*eta*gamma*sqrt( sigma*(1 + kappa) / kappa )+...
            %                  2*tau*sigma*gamma^2*( 1 + kappa ) / kappa ) / log( 1 / kappa ); 
            %
            % where:
            %
            % kappa < 1
            %
            % || u_opt ||_2 < gamma * || x_0 ||_Q
            %
            % || V(x+_0) - V(x_1) || < eta_bar * || x+_0 - x_1 ||_2 + tau_bar * || x+_0 - x_1 ||^2_2
            % eta = eta_bar * SQRT( max_eigenvalue( BETA^T * B^T * B * BETA ) )
            % tau = 0.5 * tau_bar * max_eigenvalue( BETA^T * B^T * B *BETA ) 
            % BETA = [ I, 0, ..., 0]
            
            
            % CDC 2023 formulation:
            % kappa < 1
            %
            % || u_opt ||_2 < gamma * || x_0 ||_Q
            %
            % || V(x+_0) - V(x_1) || < eta_bar_1 * || x+_0 - x_1 ||_2 + eta_bar_2 * || x+_0 - x_1 ||^2_2
            % eta_1 = eta_bar_1 * SQRT( max_eigenvalue( BETA^T * B^T * B * BETA ) )
            % eta_2 = 0.5 * eta_bar_2 * max_eigenvalue( BETA^T * B^T * B *BETA ) 
            % BETA = [ I, 0, ..., 0]
             
            if isempty(varargin)
                obj.maxiter = 5; % Default value
            else
                % Input parser
                p = inputParser;
                addOptional(p,'kappa',-Inf);
                addOptional(p,'gamma',-Inf);
                addOptional(p,'sigma',-Inf);
                addOptional(p,'eta',-Inf);
                addOptional(p,'tau',-Inf);
                parse(p,varargin{:});
                % Extract variables
                kappa   = p.Results.kappa;
                gamma   = p.Results.gamma;
                sigma   = p.Results.sigma;
                eta     = p.Results.eta;
                tau     = p.Results.tau;
                % evaluation of MAXITER
                obj.maxiter = 2*log( 2*eta*gamma*sqrt( sigma*(1 + kappa) / kappa )+...
                              2*tau*sigma*gamma^2*( 1 + kappa ) / kappa ) / log( 1 / kappa );
            end
        end

        function obj = getLargeQP(obj)
            % Get the QP paramters for large QP problems
            % This is for other algorithms to use
            obj.dmin = obj.dmin / obj.cons_mul; % goes back to true bound
            obj.dmax = obj.dmax / obj.cons_mul; % goes back to true bound
            obj.QP_H = 2*blkdiag(obj.R, kron(speye(obj.N-1),blkdiag(obj.Q,obj.R)), obj.P);
            qxu = [obj.Q*obj.xr;obj.R*obj.ur];
            obj.QP_g = -2*[obj.R*obj.ur; kron(ones(obj.N-1,1),qxu); obj.P*obj.xNr];
            % construction of equality contraint ------------------------------
            nx = obj.nx; nu = obj.nu; N = obj.N;
            KKT_con = zeros(nx*N,(nx+nu)*N);
            KKT_con(1:nx,1:nu+(nx+nu)) = [obj.H0,-obj.Gk]; % TODO: might be wrong
            for kk = 1:N-2
                KKT_con((nx)*kk+1:(nx)*(kk+1), nu+(nx+nu)*(kk-1)+1:nu+(nx+nu)*(kk+1)  ) = [obj.Hk,-obj.Gk];
            end
            kk = N-1;
            KKT_con((nx)*kk+1:(nx)*(kk+1), nu+(nx+nu)*(kk-1)+1:nu+(nx+nu)*kk + nx ) = [obj.Hk,-obj.GN];
            Aeq = KKT_con;
            leq = zeros((obj.nx)*(obj.N),1);
            ueq = zeros((obj.nx)*(obj.N),1);
            % construction of inequality constraint --------------------------------
            % Aineq = speye(obj.N *(obj.nx + obj.nu));
            Aineq_small = [obj.C, obj.D ; zeros(obj.nu,obj.nx), eye(obj.nu)];
            Aineq = blkdiag(eye(obj.nu),kron(speye(obj.N-1),Aineq_small), obj.C ); % TODO
            lineq = [obj.umin; repmat([obj.dmin;obj.umin], obj.N-1, 1); obj.dmin];
            uineq = [obj.umax; repmat([obj.dmax;obj.umax], obj.N-1, 1); obj.dmax];
            obj.QP_A = [Aeq; Aineq];
            obj.QP_l = [leq; lineq];
            obj.QP_u = [ueq; uineq];

            obj.dmin = obj.dmin * obj.cons_mul; % goes back to shrinked bound
            obj.dmax = obj.dmax * obj.cons_mul; % goes back to shrinked bound
        end

        %function obj = preCQP(obj)
        %    % pre-compute the KKT matrix
        %    % min 1/2 x'Qx + q'x; s.t. Ax=b | lam and the corresponding KKT matrix is given by
        %    % [Q,A'; * [x;   = [-q;
        %    %  A,0]     lam]    b]
        %    % The KKT matrix is given by
        %    % KKT = [Q,A';Q,0]
        %    % KKT_obj = Q
        %    tic;
        %    R = obj.R; Q = obj.C'*obj.Q*obj.C; P = obj.P;
        %    nx = obj.nx; nu = obj.nu; N=obj.N;
        %    KKT_obj = 2*blkdiag(blkdiag(R),kron(eye(N-1),blkdiag(Q,R)),P);
        %    KKT_con = zeros(nx*N,(nx+nu)*N);
        %    KKT_con(1:nx,1:nu+(nx+nu)) = [obj.H0,-obj.Gk]; % TODO: might be wrong
        %    for kk = 1:N-2
        %        KKT_con((nx)*kk+1:(nx)*(kk+1), nu+(nx+nu)*(kk-1)+1:nu+(nx+nu)*(kk+1)  ) = [obj.Hk,-obj.Gk];
        %    end
        %    kk = N-1;
        %    KKT_con((nx)*kk+1:(nx)*(kk+1), nu+(nx+nu)*(kk-1)+1:nu+(nx+nu)*kk + nx ) = [obj.Hk,-obj.GN];

        %    obj.KKT_con = KKT_con;
        %    obj.KKT = [KKT_obj KKT_con';KKT_con zeros(size(KKT_con,1))];
        %    obj.KKT_obj = sparse(KKT_obj);
        %    obj.KKT = sparse(obj.KKT);
        %    obj.dKKT = decomposition(obj.KKT);
        %    fprintf("precomputed KKT matrix for coupled QP with %f seconds\n",toc);
        %end

        %function obj = solveCoupleUpdate(obj,x0)
        %    % Solve the coupled QP problem of the following form
        %    % and update the primal and dual variable
        %    % min 1/2 x'Qx + q'x; s.t. Ax=b | lam and the corresponding KKT matrix is given by
        %    % [Q,A'; * [x;   = [-q;
        %    %  A,0]     lam]    b]
        %    ny = obj.nu+obj.nx;
        %    xi_v = vertcat(obj.xi{:}); % the vectorized xi
        %    z_v = vertcat(obj.z{:});
        %    [del_lam,primal] = pempc_solve_couple_mex(xi_v,z_v,x0);
        %    %tmp = 2*xi_v - y_v; %+ obj.yr; % TODO
        %    %obj.KKT_res = [obj.KKT_obj*tmp;obj.KKT_con_res];
        %    %sol_cQP = obj.dKKT\obj.KKT_res;
        %    %primal = sol_cQP(1:ny*obj.N);
        %    %del_lam  = sol_cQP(ny*obj.N+1:end);
        %    %obj.z{1} = primal(1:obj.nu);
        %    obj.z{obj.N+1} = primal(end-obj.nx+1:end);
        %    for kk = 1:obj.N
        %        istart = (kk-1)*(ny)+1;
        %        obj.z{kk} = primal(istart:istart+ny-1);
        %    end
        %    nxz = obj.nx;
        %    for kk = 1:obj.N+1
        %        istart = (kk-1)*(nxz)+1;
        %        obj.del_lam{kk} = del_lam(istart:istart+nxz-1);
        %        obj.lam{kk} = obj.lam{kk} + obj.del_lam{kk};
        %    end
        %end

    end % end of public method
end

