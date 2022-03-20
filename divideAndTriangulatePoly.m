function [polys,bad] = divideAndTriangulatePoly(parent)
%divideAndTriangulatePoly will fracture a polygon into 3 or 4 sided shapes.
%
%[polys] = divideAndTriangulatePoly(parent) divides a parent polygon (a
%struct with fields x and y that are column vectors) and then triangulates
%the resulting polygons if needed to have all 3 or 4-sided shapes.

%If we have a very complicated polygon, try to split it with
%DIVIDEXY, so we avoid too many nasty triangles.
if length(parent.x>10)
    ndiv=round(length(parent.x)/4);
    polys=DIVIDEXY(parent,ndiv,ndiv);
    polys=polys(:);
    [polys,bad]=checkPolys(polys,parent);
    polys=fixPolys(polys);
    
    %Go through the resulting polygons and triangulate if needed.
    polys2={};
    for ipp=1:length(polys)
        if length(polys{ipp}.x)>4
            pp=triangulatePoly(polys{ipp});
            [pp,bad]=checkPolys(pp,polys{ipp});
            pp=fixPolys(pp);
            [pp,bad]=checkPolys(pp,polys{ipp});
            polys2=[polys2 pp];
        else
            polys2=[polys2 polys(ipp)];
        end
    end
    polys=polys2;
else %not so many vertices, so just triangulate
    polys=triangulatePoly(parent);
    [polys,bad]=checkPolys(polys,parent);
    polys=fixPolys(polys);
    [polys,bad]=checkPolys(polys,parent);

end


end

