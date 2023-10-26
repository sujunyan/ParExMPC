function obj = preRiccati(obj)
% The Factorization phase of the Riccati method
obj.Ric_P{obj.N+1} = 2*obj.P;
for n = obj.N:-1:1
    AtPA = obj.A'*obj.Ric_P{n+1}*obj.A;
    BtPB = obj.B'*obj.Ric_P{n+1}*obj.B;
    BtPA = obj.B'*obj.Ric_P{n+1}*obj.A;
    obj.Ric_Lam{n} = chol(2*obj.R+BtPB,'lower');
    obj.Ric_Lam_inv{n} = inv(obj.Ric_Lam{n});
    obj.Ric_L{n} = obj.Ric_Lam{n} \ BtPA;
    PP = 2*obj.Q + AtPA - obj.Ric_L{n}'* obj.Ric_L{n};
    obj.Ric_P{n} = 1/2 * (PP+PP');
end
end