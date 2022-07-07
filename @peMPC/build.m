function obj = build(obj)
    obj = obj.init; % initialize the pempc solver
    obj = obj.getMPTfunc; % generate mpt codes
    obj.toC; % generate c code
    obj.compile_c;
end