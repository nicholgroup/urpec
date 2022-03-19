function [polys] = fixPolys(polys)
%fixPolys calls fixPoly for a cell array of polygons.

for ip=1:length(polys)
    [polys{ip}.x polys{ip}.y]=fixPoly(polys{ip}.x,polys{ip}.y);
end
end

