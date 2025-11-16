classdef RRLBMPC
    % parallel explicit MPC controller
    % TODO: refer to the paper: xxxxxxx
    % TODO: add more information on the RRLB related properties
    % Author:         Junyan Su
    % The pempc object to solve the problem of the form
    %   \min \sum (x[k]-x_r)' Q (x[k]-x_r) for k = 1...N-1
    %           + (u[k]-u_r)' R (u[k]-u_r)
    %           + (x[N]-x_Nr)' P (x[N]-x_Nr)
    %   s.t. x[k+1] = Ax[k] + Bu[k]
    %        Cx * x[k] <= dx
    %        Cu * u[k] <= du
    % Optional input:
    %   xr: the reference state. The default value is zero
    %   xNr: the reference terminal state. The default value is zero
    %   ur: the reference control input. The default value is zero
    %   xmin,xmax: the state constraint. The default value is [-inf,inf]
    %   umin,umax: the control input constraint. The default value is [-inf,inf]
    %   N: the time horizon. The default value is 10
    %   cons_mul: the constraint multiplier, often set less than 1 to avoid constraint violation. The default value is 1
    %   par_flag: the boolean flag to enable/disable the parallel computing. The default value is true
    %   par_threshold: enable the parallel computing if the time horizon is larger than the threshold. The default value is 20
    % Call:
    %       mpc0 = peMPC(A,B,Q,R,P,{,optional inputs})

    properties
        % The System dynamics
        % x[k+1] = A x[k] + B u[k]
        % Cx * x[k] <= dx
        % Cu * u[k] <= du
        A
        B
        Cx
        dx
        Cu
        du
        % TODO: for now, we assume the system is fully observable

        % Reference values
        % TODO: for now, we treat it as constant
        xr          % reference state trajectory
        ur          % reference input trajectory
        xNr         % reference terminal state
        zr          % The reference stack variable, have size (nx+nu)*N

        % Objective function
        Q           %  (x[k]-x_r[k])' Q (x[k]-x_r[k]) for k = 1...N-1
        R           %  (u[k]-u_r[k])' R (u[k]-u_r[k])
        P           %  (x[N]-x_r[N])' P (x[N]-x_r[N])
        Sigma0      % blkdiag(R,S)
        Sigmak      % blkdiag(Q,R,S)
        NSigmak     % kron(eye(N-1),Sigmak)
        SigmaN      % P

        % RRLB related 
        wx          % weight for RRLB on state. 
        delta_relax     % The tolerance in the relaxed log barrier. If x < delta_relax, then the function becomes a quadratic function.
        rho         % The weight for the relaxed log barrier function. The final stage cost is L(x,u) = l(x,u) + rho * bx(x)


        % the dimentions
        nx
        nu
        mx  % number of state constraints
        mu  % number of input constraints

        % controller parameters
        maxiter     % the maximum iteration
        N           % the time horizon
        isFirst     % the flag to indicate that if this object has ben used. For the first time, we run a large number of itertions in aladin to initilize the MPC controller.

        % the stored variable... For the new interface, we should not store any internal variables...
        % lam         % the Lagrangian multiplier
        % del_lam     % TODO: for testing
        % xi          % the alternating direction

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

        % the stack variables in the condensed form
        compactA
        compactb
        compactQ
        compactR

        % z1          % the stacked variable for states
        % z2          % the stacked variable for control inputs

        % ALADIN related properties -----------------------
        ALADIN_H1
        ALADIN_H2

        % ADMM related properties (for comparison only) ------------------------------
        ADMM_sigma   % the penalty parameter in ADMM



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Below are deprecated properties %
        %%%
        
        % State and input constraints
        % dmin        % dmin <= C*x + Du <= dmax
        % dmax
        % umin        % umin <= u[k] <= umax
        % umax
        % z           % The stacked variable [x;u]
        % zmin
        % zmax

        % The converted system dynamics---in case we have a different system format
        % G[k+1]z[k+1] = H[k]z[k] + h[k]
        % where z[0] = u[0]; z[k] = [x[k]; u[k]]; z[N] = x[N]
        % GN
        % Gk
        % NGkt        % kron(eye(N-1),Gk')
        % GN
        % Hk
        % NHkt        % kron(eye(N-1),Hk')
        % H0

    end % End of the properties

    methods (Access = public)

        function out = isPosDef(obj,A)
            % determine if matrix A is positive definite.
            tol = 1e-8;
            d = eig(A);
            out = all(d > tol);
        end

        function obj = RRLBMPC(A,B,Q,R,P,varargin)
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
            %   Cx, dx: the constraint Cx * x[k] <= dx
            %   Cu, du: the constraint Cu * u[k] <= du
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

            addOptional(p,'Cx', []); addOptional(p,'dx',[]);
            addOptional(p,'Cu', []); addOptional(p,'du',[]);

            addOptional(p,'rho', 1e-4);


            addOptional(p,'N',10);
            addOptional(p,'cons_mul',1);
            addOptional(p,'mptSolver','plcp');
            addOptional(p,'par_flag', true);
            addOptional(p,'par_threshold', 20);

            addOptional(p,'delta', 1e-2);
            addOptional(p,'ADMM_sigma', 1.0);

            parse(p,varargin{:});
            obj.xr = p.Results.xr;
            obj.xNr = p.Results.xNr;
            obj.ur = p.Results.ur;

            obj.Cx = p.Results.Cx; obj.dx = p.Results.dx;
            obj.mx = size(obj.Cx,1);
            obj.Cu = p.Results.Cu; obj.du = p.Results.du;
            obj.mu = size(obj.Cu,1);
           
            obj.ADMM_sigma = p.Results.ADMM_sigma;
            obj.delta_relax = p.Results.delta;

            obj.N = p.Results.N;
            obj.cons_mul = p.Results.cons_mul;
            obj.rho = p.Results.rho;

            % obj.dmin = obj.dmin * obj.cons_mul;
            % obj.dmax = obj.dmax * obj.cons_mul;
            obj.mptSolver = p.Results.mptSolver;
            obj.use_parallel = p.Results.par_flag;
            obj.parallel_threshold = p.Results.par_threshold;


        end

        function obj = getWx(obj)
            % Get the weight vector for RRLB
            % We want to solve the following simple QP
            % min_w  ||w - 1||_2^2
            % s.t. \sum_r w_r * a_r == 0
            % here a_r is the derivative of the r-th constraint at the point zero.

            a_vec = zeros(obj.mx, 1);
            delta = obj.delta_relax;
            for r = 1:obj.mx
                % get the derivative of the r-th constraint at zero
                dr_i = obj.dx(r);
                if dr_i > delta
                    ai = - 1 / dr_i;
                else
                    ai = (dr_i - 2*delta) / delta;
                end
                a_vec(r) = ai;
            end

            % By solving the KKT condition, we have
            lam_tmp = 2 * sum(a_vec) / (a_vec' * a_vec);
            
            obj.wx = 1 - 0.5 * lam_tmp * a_vec;



        end

        function obj = getCompactA(obj)
            %%% Create the compact A matrix for the condensed form
            % Generated by Gemini, need to verify
            A = obj.A; B = obj.B; N = obj.N;


            [m, ~] = size(A);
            [~, p] = size(B);
            % Pre-allocate cell array for blocks
            block_cells = cell(N, N);
    
            % Store powers of A * B for reuse
            AB_powers = cell(N, 1);
            current_A_power_B = B; % A^0 * B
            AB_powers{1} = current_A_power_B;
            for k = 2:N
                current_A_power_B = A * current_A_power_B; % Calculate A^(k-1) * B
                AB_powers{k} = current_A_power_B;
            end
            % Fill the blocks
            for j = 1:N % Block columns
                for i = j:N % Block rows (lower triangular)
                    block_cells{i, j} = AB_powers{i - j + 1};
                end
                % Fill upper triangular part with zeros (represented by empty cells for cell2mat or explicit zero matrices)
                for i = 1:(j-1)
                    block_cells{i, j} = zeros(m, p);
                end
            end
    
            % Convert cell array of blocks to a single matrix
            obj.compactA = cell2mat(block_cells);

        end


        function obj = getCompactForm(obj)
            % Get the compact form for the problem
            % This is mainly used for verification
            % Minimize f1(z1) + f2(z2)
            %   s.t. z1 = A z2 + bb

            
            obj = obj.getCompactA;

            obj.compactb = zeros(obj.N * obj.nx, 1);

            obj.compactQ = blkdiag(kron(speye(obj.N-1), obj.Q), obj.P);
            obj.compactR = kron(speye(obj.N), obj.R);
            
        end

        function obj = updateX0(obj, x0)
            %%% Update according to the current state x0

            n = obj.nx;
            current_A_power = eye(n); % Represents A^0
            for k = 1:obj.N
                current_A_power = current_A_power * obj.A; % Calculate A^k
                block_k = current_A_power * x0;        % Calculate A^k * x0
                obj.compactb(((k-1)*n + 1) : (k*n), 1) = block_k; % Place block into b
            end
        end

        function res = getZ1k(obj, z1, k)
            % get the z1_k from the stacked variable z1
            res = z1((k-1)*(obj.nx)+1 : (k-1)*(obj.nx)+obj.nx);
        end

        function res = getZ2k(obj, z2, k)
            % get the z2_k from the stacked variable z2
            res = z2((k-1)*(obj.nu)+1 : (k-1)*(obj.nu)+obj.nu);
        end
      
        function res = f1(obj, z1)
            res = z1' * obj.compactQ * z1;
            delta = obj.delta_relax;
            for k = 1:obj.N-1
                z1_k = getZ1k(obj, z1, k);

                rrlb = RRLB(obj, delta, obj.rho, obj.wx, obj.Cx, obj.dx, z1_k);
                res = res + rrlb;
            end

        end


        function res = f2_no_cons(obj, z2)
            res = z2' * obj.compactR * z2;
        end
        

        % Other help methods go here ----------------------------------------

        function res = relax_barrier(obj, delta, x)
            % The relaxed log barrier function
            res = 0;
            if x >= delta
                res = -log(x);
            else
                res = 0.5 * (( (x - 2*delta)/delta)^2 - 1) - log(delta);
            end
        end

        function res = relax_barrier_d2(obj, delta, x)
            % The second derivative of the relaxed log barrier function
            res = 0;
            if x >= delta
                res = 1 / (x^2);
            else
                res = 1 / (delta^2);
            end
        end

        function res = RRLB(obj, delta, rho, wx, Cx, dx, x)
            % The relaxed recentered log barrier (RRLB) function
            res = 0;
            for r = 1:obj.mx
                rlb1 = relax_barrier(obj, delta, dx(r) - Cx(r,:) * x);
                rlb0 = relax_barrier(obj, delta, dx(r));
                res = res + rho * wx(r) * (rlb1 - rlb0);
            end
        end


    end % end of public method
end

