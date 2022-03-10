function [bad] = checkArea(poly,polys2Add)
%CHECKAREA checks that the area of a polygon is equal to the area of the fractured polygons
%
% CHECKAREA computes the area of a polygon (poly) and adds the areas of the
% fractured polygons (polys2Add). The output argument (bad) is 1 if the
% areas are not the same. This can happen if the polygon is fractured into
% shapes that are not good polygons.

bad=0;
a1=polyarea(poly.x,poly.y);
a2=0;
for j=(1:1:length(polys2Add))
    a2=a2+polyarea(polys2Add{j}.x,polys2Add{j}.y);
end

if ~(round(a1,6)==round(a2,6))
    bad=1;
end
end

