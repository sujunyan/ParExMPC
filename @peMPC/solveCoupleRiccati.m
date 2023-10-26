function obj = solveCoupleRiccati(obj,x0)
% Solve the coupled QP problem of the following form with Riccati backward-farward sweep strategy.
% Refer to the Alg.1-2 of the paper:
% Efficient Implementation of the Riccati Recursion for Solving Linear-Quadratic Control Problems
tic;
p{obj.N+1} = - 2*obj.P * (2*obj.xi{obj.N+1}-obj.z{obj.N+1});
for n = obj.N:-1:1
    tmp = (2*obj.xi{n}-obj.z{n});
    qn = - 2*obj.Q * tmp(1:obj.nx);
    sn = - 2*obj.R * tmp(obj.nx+1:end);
    l{n} = obj.Ric_Lam_inv{n} * (sn + obj.B'*p{n+1});
    p{n} = qn + obj.A'*p{n+1} - obj.Ric_L{n}'*l{n};
end
xn = x0;
xix0 = obj.xi{1};
%tmp_lam{1} = -obj.Q*(3*x0-4*xix0(1:obj.nx));
tmp_lam{1} = obj.Ric_P{1}*xn+p{1};
for n = 1:obj.N
    un = (- obj.Ric_Lam_inv{n}')*(obj.Ric_L{n} * xn+l{n});
    tmp_z{n} = [xn;un];
    xn = obj.A*xn + obj.B*un;
    tmp_lam{n+1} = obj.Ric_P{n+1}'*xn + p{n+1};
    %fprintf("xn:%f,%f\n",xn);
end
tmp_z{obj.N+1} = xn;
obj.z = tmp_z;
for ii = 1:obj.N+1
    obj.lam{ii} = obj.lam{ii} + tmp_lam{ii};
end
end