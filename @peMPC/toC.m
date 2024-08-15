function toC(obj)
% peMPC method that generate a branch of C files for one sampling time of MPC
class_folder =  fileparts(mfilename('fullpath'));
c_folder = [class_folder, filesep,'c_code'];
%copyModifyMPTc;
writeHeadFile(obj,c_folder);
end

function writeHeadFile(obj,folder)
% write ther headers that stores the pre-computed constant variables
fileID = fopen([folder,filesep,'MPC_problem.h'],'w+');
fprintf(fileID,"//Generated header that stores the pre-defined varibales\n");
fprintf(fileID,"#ifndef __PREDEFINED_H__\n#define __PREDEFINED_H__\n");
fprintf(fileID,"#include <stddef.h>\n");
fprintf(fileID,"static const size_t nx = %d;\n",obj.nx);
fprintf(fileID,"static const size_t nu = %d;\n",obj.nu);
fprintf(fileID,"static const size_t N = %d;\n",obj.N);
fprintf(fileID,"#define PEMPC_N %d // for control the openmp\n", obj.N);
if obj.use_parallel
    fprintf(fileID,"#define USE_OPENMP (PEMPC_N > %d) // for control the openmp\n", obj.parallel_threshold);
else
    fprintf(fileID,"#define USE_OPENMP 0 // for control the openmp\n");
end
printMatrixAsCstring(obj.Q,'Q','fileID',fileID);
printMatrixAsCstring(obj.P,'P','fileID',fileID);
printMatrixAsCstring(obj.R,'R','fileID',fileID);
printMatrixAsCstring(obj.A,'A','fileID',fileID);
printMatrixAsCstring(obj.B,'B','fileID',fileID);
printMatrixAsCstring(obj.ur,'ur','fileID',fileID);
printMatrixAsCstring(obj.xr,'xr','fileID',fileID);
printMatrixAsCstring(obj.xNr,'xNr','fileID',fileID);
printCellAsCstring(obj.Ric_Lam,'Ric_Lam','fileID',fileID);
printCellAsCstring(obj.Ric_Lam_inv,'Ric_Lam_inv','fileID',fileID);
printCellAsCstring(obj.Ric_L,'Ric_L','fileID',fileID);
printCellAsCstring(obj.Ric_P,'Ric_P','fileID',fileID);
fprintf(fileID,"#endif\n");
fclose(fileID);
end


function printCellAsCstring(Ce,name,varargin)
% print to fileID as array in c-type
% column-major
% return ----------------------------
% a cell of strings where each cell
p = inputParser;
addOptional(p,'precision','double');
addOptional(p,'fileID',1);
parse(p,varargin{:});
precision = p.Results.precision;
fileID = p.Results.fileID;

fprintf(fileID,['static ', precision,' ',name, '[] = {\n']);
for ii = 1:length(Ce)
    A = Ce{ii};
    [nRow,nCol] = size(A);
    for col = 1:nCol
        str = [];
        for row = 1:nRow
            if (row == nRow) && (col == nCol) && (ii==length(Ce))
                str = [str, sprintf(' %.15e', A(row,col))];
            else
                str = [str, sprintf(' %.15e,\t', A(row,col))];
            end
        end
        fprintf(fileID,[str,'\n']);
    end
end

fprintf(fileID,'};\n');
end

function s = printMatrixAsCstring(A,name,varargin)
% print to fileID as array in c-type
% column-major
% return ----------------------------
% a cell of strings where each cell
p = inputParser;
addOptional(p,'precision','double');
addOptional(p,'fileID',1);
parse(p,varargin{:});
precision = p.Results.precision;
fileID = p.Results.fileID;

[nRow,nCol] = size(A);
s{1} = ['static ',precision,' ',name, '[] = {'];
for col = 1:nCol
    str = [];
    for row = 1:nRow
        if (row == nRow) && (col == nCol)
            str = [str, sprintf(' %.15e', A(row,col))];
        else
            str = [str, sprintf(' %.15e,\t', A(row,col))];
        end
    end
    s{col+1} = str;
end
s{nCol+2} = '};';
for ii = 1:nCol+2
    fprintf(fileID,[s{ii},'\n']);
end
end