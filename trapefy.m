function [polys] = trapefy(poly)
%TRAPEFY fractures a polygon into trapezoids.
%   
%   polys=TRAPEFY(polys) is a cell array of polygons, where each of the
%   polygons is a trapezoid, and the union of polys is the original polygon.
%   Poly is a struct with fields x and y, which specify the vertices of hte
%   polygon.

poly=polyshape(poly.x,poly.y);

%The y coordinates of all of the vertices.
v=poly.Vertices(:,2);

polys=poly;

%Split the polygon horitonzally at all of the y coordinates found above.
for iv=1:size(poly.Vertices,1)   
        
    for ip=length(polys):-1:1       
        p=struct();
        p.x=polys(ip).Vertices(:,1);
        p.y=polys(ip).Vertices(:,2);
        
        [p1 p2]=DIVIDEY(p,v(iv));

        if ~isempty(p1)
            poly1=polyshape(p1.x,p1.y);
            polys1=regions(poly1);
        else
            polys1=[];
        end     
        if ~isempty(p2)
            poly2=polyshape(p2.x,p2.y);
            polys2=regions(poly2);
        else
            polys2=[];
        end        
        polys=[polys(1:ip-1); polys1; polys2; polys(ip+1:end)];
    end
end

%convert it to a cell array with the right format.
pc=polys;
polys={};
for ip=1:length(pc)
    polys{ip}.x=pc(ip).Vertices(:,1);
    polys{ip}.y=pc(ip).Vertices(:,2);
end

end



