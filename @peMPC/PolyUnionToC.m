function PolyUnionToC(obj, fid, PU, prefix)
% Modified from original method @PolyUnion.toC() form MPT Toolbox.
% We modify this file specific for our application for better performance.
%
% Exports PWA/PWQ function to a C-code that uses a sequential search to
% evaluate the function.
%
%   syntax: U.toC('function', 'output','tie_break_fcn')
%
%   'fid' - fileID to wirte to
%   'PU' - the PolyUnion object to export
%   'prefix' - the prefix to identify different data
%
%  Generated code contains a single .h file which contains the necessary matrix
%  for evaluation. We implement the evaluation function inside the pempc.c file.
%  Note that we have 3 PWA functions in total to output, so in this file, we
%  only produces the data and leave other macro to other function

global MPTOPTIONS

PU = PU{1};
function_name = 'primal';
tie_break_fcn = 'obj';

% check that all PolyUnions have the same list of functions ----------------------
lF = PU(1).listFunctions;
for i=2:numel(PU)
    for j=1:numel(lF)
        if ~PU(i).hasFunction(lF{j})
            error('All objects in the array must have attached the same list of functions.');
        end
    end
end

% check that the PolyUnion is full-dimensional
for i=1:numel(PU)
    if ~PU(i).isFullDim
        error('The export to C-code works only for full-dimensional partitions.');
    end
end

% only PWA/PWQ functions are allowed
for i=1:numel(PU)
    % is the requested function present?
    if ~PU(i).hasFunction(function_name)
        error('No such function "%s" in the object.', function_name);
    end
    for j=1:PU(i).Num
        if ~(isa(PU(i).Set(j).Functions(function_name), 'AffFunction') || ...
                isa(PU(i).Set(j).Functions(function_name), 'QuadFunction') )
            error('Only quadratic and affine functions are supported.');
        end
        if ~isempty(tie_break_fcn)
            if ~(isa(PU(i).Set(j).Functions(tie_break_fcn), 'AffFunction') || ...
                    isa(PU(i).Set(j).Functions(tie_break_fcn), 'QuadFunction') )
                error('Only quadratic and affine tie-break functions are supported.');
            end
        end
    end
end

% check that all PolyUnions have requested functions and the range is the same
for i=1:numel(PU)
    R = PU(1).Set(1).Functions(function_name).R;
    for j=1:PU(i).Num
        if R~=PU(i).Set(j).Functions(function_name).R
            error('The requested function "%s" must have the same range in all objects in the array.',function_name);
        end
    end
end

% single or double precision to export?
precision = MPTOPTIONS.modules.geometry.unions.PolyUnion.toC.precision;
precision = strtrim(lower(precision));
if isempty(precision)
    precision = 'double';
elseif ~isequal(precision,'single') && ~isequal(precision,'double')
    error('The specified precision in the option can be either "single" or "double".');
end
if isequal(precision,'single')
    precision = 'float';
end

% extract polyhedra with control law
Pn = [PU.Set];
nr = [PU.Num];
total_nr = sum(nr);

% extract hyperplane representation
An = cell(total_nr,1);
bn = cell(total_nr,1);
[An{:}]=deal(Pn.A);
[bn{:}]=deal(Pn.b);
if ~iscell(An),
    An = {An};
    bn = {bn};
end

% count number of constraints
nctotal = 0;
for ii=1:total_nr,
    nctotal = nctotal + size(Pn(ii).H,1);
end

% extract dimensions
nx = PU(1).Dim; % domain
nu = PU(1).Set(1).Functions(function_name).R; % range

Fi = cell(total_nr,1);
Gi = Fi;
for i=1:total_nr
    Fi{i}=Pn(i).Functions(function_name).F;
    Gi{i}=Pn(i).Functions(function_name).g;
end

% extract tie-break function
HTBi = cell(total_nr,1);
FTBi = cell(total_nr,1);
GTBi = FTBi;
for i=1:total_nr
    HTBi{i}=Pn(i).Functions(tie_break_fcn).H;
    FTBi{i}=Pn(i).Functions(tie_break_fcn).F;
    GTBi{i}=Pn(i).Functions(tie_break_fcn).g;
end

%% write mpt_getInput.c --------------------------------------------

fprintf(fid,'/* Generated on %s by MPT %s */ \n\n',datestr(now), MPTOPTIONS.version);

if ~isempty(tie_break_fcn)
    fprintf(fid, '#include <float.h>\n\n');
end
fprintf(fid, '#define %sMPT_NR (%d)\n', prefix, total_nr);

% write inequality constraints A*x <= b for each polytope
ctr = 0;
fprintf(fid, '\nstatic %s %sMPT_A[] = {\n', precision,prefix);
writeCellMatrices(An,fid,precision);
%str_buffer = "";
%for ii = 1:total_nr,
%    Ai = An{ii};
%    nc = size(Ai, 1);
%    for jj = 1:nc,
%        a = Ai(jj, :);
%        for kk = 1:length(a),
%            ctr = ctr + 1;
%            if ctr<nctotal*nx,
%                if isequal(precision,'float')
%                    str_buffer = str_buffer + sprintf('%.7e,\t',a(kk));
%                    %fprintf(fid, '%.7e,\t', a(kk));
%                else
%                    str_buffer = str_buffer + sprintf('%.14e,\t',a(kk));
%                    %fprintf(fid, '%.14e,\t', a(kk));
%                end
%            else
%                if isequal(precision,'float')
%                    str_buffer = str_buffer + sprintf('%.7e ',a(kk));
%                    %fprintf(fid, '%.7e ', a(kk));
%                else
%                    str_buffer = str_buffer + sprintf('%.14e ',a(kk));
%                    %fprintf(fid, '%.14e ', a(kk));
%                end
%            end
%            if mod(ctr, 5)==0,
%              %str_buffer = str_buffer + '\n';
%            end
%        end
%    end
%end
%str_buffer = str_buffer + '};\n\n';

%ctr = 0;
fprintf(fid, 'static %s %sMPT_B[] = {\n',precision,prefix);
writeCellMatrices(bn,fid,precision);
%for ii = 1:total_nr,
%    bi = bn{ii};
%    nc = size(bi, 1);
%    for jj = 1:nc,
%        ctr = ctr + 1;
%        if ctr<nctotal,
%            if isequal(precision,'float')
%                fprintf(fid, '%.7e,\t', bi(jj));
%            else
%                fprintf(fid, '%.14e,\t', bi(jj));
%            end
%        else
%            if isequal(precision,'float')
%                fprintf(fid, '%.7e ', bi(jj));
%            else
%                fprintf(fid, '%.14e ', bi(jj));
%            end
%        end
%        if mod(ctr, 5)==0,
%            fprintf(fid, '\n');
%        end
%    end
%end
%fprintf(fid, '};\n\n');

fprintf(fid, 'static int %sMPT_NC[] = {\n',prefix);
for ii = 1:total_nr,
    if ii < total_nr,
        fprintf(fid, '%d,\t', size(Pn(ii).H,1));
    else
        fprintf(fid, '%d ', size(Pn(ii).H,1));
    end
    if mod(ii, 5)==0,
        fprintf(fid, '\n');
    end
end
fprintf(fid, '};\n\n');

F_name = [prefix, 'MPT_F'];
G_name = [prefix, 'MPT_G'];
sub_write_matrix(Fi, F_name, fid, precision);
sub_write_matrix(Gi, G_name, fid, precision);
sub_write_matrix(HTBi, [prefix, 'MPT_HTB'], fid, precision);
sub_write_matrix(FTBi, [prefix, 'MPT_FTB'], fid, precision);
sub_write_matrix(GTBi, [prefix, 'MPT_GTB'], fid, precision);

end % end of PolyUnionToC

function sub_write_matrix(matrices, name, fid, precision)
%
% writes a cell of matrices to a file with a given name
%
% syntax:
%         matrix - a cell array with matrix data
%         name - name of the matrix given as string
%         fid - file identificator
%         precision - either "float" or "double"

nr = numel(matrices);
[nu,nx] = size(matrices{1});

nctotalh = nu*nx*nr;
ctr = 0;
%fprintf(fid, '\n');
fprintf(fid, '\nstatic %s %s[] = {\n', precision, name);
writeCellMatrices(matrices,fid,precision);
%An_mat = cell2mat(matrices);
%An_v = reshape(An_mat',[],1);
%buf2 = sprintf("%.14e,\t",An_v(1:end-1));
%buf2 = buf2 + sprintf("%.14e ",An_v(end)) + '};\n\n';
%fprintf(fid,buf2);
%for ii = 1:nr,
%    M = matrices{ii};
%    for jj = 1:nu,
%        h = M(jj, :);
%        for kk = 1:nx,
%            ctr = ctr + 1;
%            if ctr<nctotalh,
%                if isequal(precision,'float')
%                    str_buffer = str_buffer +  sprintf("%.7e,\t", h(kk));
%                    %fprintf(fid, '%.7e,\t', h(kk));
%                else
%                    str_buffer = str_buffer +  sprintf("%.14e,\t", h(kk));
%                    %fprintf(fid, '%.14e,\t', h(kk));
%                end
%            else
%                if isequal(precision,'float')
%                    str_buffer = str_buffer +  sprintf("%.14e ", h(kk));
%                    %fprintf(fid, '%.7e ', h(kk));
%                else
%                    str_buffer = str_buffer +  sprintf("%.14e ", h(kk));
%                    %fprintf(fid, '%.14e ', h(kk));
%                end
%            end
%            if mod(ctr, 5)==0
%                %str_buffer = str_buffer + '\n';
%                %fprintf(fid, '\n');
%            end
%        end
%    end
%end
%str_buffer = str_buffer + '};\n\n';
%fprintf(fid,str_buffer);
end

function writeCellMatrices(matrices,fid,precision)
%% write the tailor of a cell to fid
% syntax:
%         matrix - a cell array with matrix data
%         name - name of the matrix given as string
%         fid - file identificator
%         precision - either "float" or "double"
An_mat = cell2mat(matrices);
An_v = reshape(An_mat',[],1);
if length(An_v) == 1 % there is no colon for single element array
    buf2 = sprintf("%.7e ",An_v(end)) + '};\n\n';
else
    if isequal(precision,'float')
        buf2 = sprintf("%.7e,\t",An_v(1:end-1));
        buf2 = buf2 + sprintf("%.7e ",An_v(end)) + '};\n\n';
    else
        buf2 = sprintf("%.14e,\t",An_v(1:end-1));
        buf2 = buf2 + sprintf("%.14e ",An_v(end)) + '};\n\n';
    end
end
fprintf(fid,buf2);
end