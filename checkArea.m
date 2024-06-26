function [bad] = checkArea(polys,parent,minArea)
%CHECKAREA checks that the area of a polygon is equal to the area of the fractured polygons
%
% CHECKAREA computes the area of a polygon (parentPoly) and adds the areas of the
% fractured polygons (polys). The output argument (bad) is 1 if the
% areas are not the same. This can happen if the polygon is fractured into
% shapes that are not good polygons.


bad=0;
a1=polyarea(parent.x,parent.y);
a2=0;
for j=(1:1:length(polys))
    a2=a2+polyarea(polys{j}.x,polys{j}.y);
end


if ~exist('minArea','var')
    minArea=1e-5; %suitable for when the units are microns
end


if abs(a1-a2)>minArea
    bad=1;
end
end

