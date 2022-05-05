function PulseExportGradientTable(pulse, filename, isExportX, isExportY, isExportZ)
% SYNTAX: 
% 
%   PulseExportGradientTable(pulse, filename, isExportX, isExportY, isExportZ)
%
% Exports the gradient shapes from the given pulse.
% Inputs:
%
% Variable Name    Description
% pulse            Input RF pulse structure
% filename         Export filename (without extension).
% isExportX,Y,Z    Either 0 or 1. If set to 1, the gradient along that
%                  axis will be exported as a file.

% Check filename extension. Remove if present. Append .h extension.
filenameNoExt = RemoveExtensionFromFilename(filename);
filenameNoExtX = [filenameNoExt,'X'];
filenameNoExtY = [filenameNoExt,'Y'];
filenameNoExtZ = [filenameNoExt,'Z'];
filenameX = [filenameNoExtX,'.h'];
filenameY = [filenameNoExtY,'.h'];
filenameZ = [filenameNoExtZ,'.h'];

if isExportX, ExportSubRoutine(pulse.Gx, filenameX); end
if isExportY, ExportSubRoutine(pulse.Gy, filenameY); end
if isExportZ, ExportSubRoutine(pulse.Gz, filenameZ); end


function ExportSubRoutine(gradientShape, filename)
% Normalize gradient shape to [-1, 1] interval 
maxGrad = max(abs(gradientShape));
if maxGrad~=0
    gradientShape = gradientShape./maxGrad;
end
varName = RemoveExtensionFromFilename(filename);

% Open file
fid = fopen(filename,'w');
for idx=1:numel(gradientShape)
    fprintf(fid,'%s[%d] = float(%.6f);\n', varName, idx-1, gradientShape(idx));
end
fclose(fid);


