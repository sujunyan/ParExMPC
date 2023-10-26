function obj = build(obj)
obj = obj.init; % initialize the peMPC solver
obj = obj.getMPTfunc; % generate MPT codes
obj.toC; % generate C code
obj.compile_c;
end