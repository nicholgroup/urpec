function [] = rotatePattern(fieldsFile,angles)
% Rotates a pattern by angle degrees
% fieldsFile: a "fields.mat" file produced by urpec
% angles: an array of angles in degrees to rotate a pattern.
%
% This function saves .dc2 and _fields.mat files for each rotated pattern.
%
% If you call the function like this: rotatePattern('',linspace(-1,1,12)),
% it will prompt you to load a file, and then it will make 12 copies, with
% rotation angles going from -1 to 1 degrees in 12 steps.

if ~exist(fieldsFile,'var')
    [fieldsFile pathname]=uigetfile('*fields.mat');
    dir=pwd;
    cd(pathname);
end

load(fieldsFile)

polygons=fields.polygons;

for a=1:length(angles)
    
    rotatedPolys=polygons;
    
    R=[cosd(angles(a)) -sind(angles(a)); sind(angles(a)) cosd(angles(a))];
    for i=1:length(rotatedPolys)
        for j=1:size(rotatedPolys(i).p,1)
            rotatedPolys(i).p(j,:)=(R*rotatedPolys(i).p(j,:)');
        end
        
    end
    
    fieldsFileName=[fieldsFile(1:end-11) strrep(sprintf('_R%1.1f',angles(a)),'.','_') '_fields.mat'];
   
%     cadName=fields.cadFile;
%     dc2FileName=[cadName(1:end-4) strrep(sprintf('_R%1.1f',angles(a)),'.','_') '.dc2'];
%     dc2write(rotatedPolys,dc2FileName);
%     fields.cadFile=dc2FileName;


    fields.polygons=rotatedPolys;
    save(fieldsFileName,'fields');
    
end

cd(dir);


end

