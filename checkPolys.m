function [polys,bad] = checkPolys(polys,parent,minArea)
%CHECKPOLYS checks fractured polygons.
%
% [polys,bad] = CHECKPOLYS(parent,polys) checks a cell array of fracture
% polygons (polys) and removes zero-area or empty entries. It also checks
% that the sum of the fractured areas equals the area of the parent
% polygon. The boolean output bad is true if the areas do not match.

if ~exist('minArea','var')
    minArea=1e-5; %suitable for when the units are microns
end

for j=(length(polys):-1:1)
    
    if isempty(polys{j}) || any(size(polys{j}.x)~=size(polys{j}.y)) || any(size(polys{j}.x)==0) 
        polys(j)=[];
    end
    
    try
        if polyarea(polys{j}.x,polys{j}.y)<minArea
            polys(j)=[];
        end
    end
    
end

%Check that the area of the fractured polygons adds up to the
%original area
bad=checkArea(polys,parent);

end

