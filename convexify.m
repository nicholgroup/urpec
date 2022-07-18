function [polys] = convexify(poly)
%CONVEXIFY fractures a polygon into convex shapes.
%   
%   polys=CONVEXIFY(polys) is a cell array of polygons, where each of the
%   polygons is convex, and the union of polys is the original polygon.
%   Poly is a struct with fields x and y, which specify the vertices of hte
%   polygon.

poly=polyshape(poly.x,poly.y);

iter=0;
convex=0;
polys=poly;

while ~all(convex)
    iter=iter+1;
    
    convex=0;
    
    for ip=length(polys):-1:1
        %Change to CW
        [x2,y2]=poly2cw(polys(ip).Vertices(:,1),polys(ip).Vertices(:,2));
        p=polyshape(x2,y2);
        
        %Find non-convex vertices
        [f,v]=checkConvex(p.Vertices(:,1)',p.Vertices(:,2)');
        convex(ip)=f;
        
        if f
            continue
        end
        
        inds=(v==1);
        i1=find(inds>0); i1=i1(1); %index of first non-convex vertex.
        
        p=struct();
        p.x=x2;
        p.y=y2;
        
        %Split along first non-convex vertex
        [p1 p2]=DIVIDEY(p,y2(i1));
        
        poly1=polyshape(p1.x,p1.y);
        poly2=polyshape(p2.x,p2.y);
        
        polys=[polys(1:ip-1); regions(poly1); regions(poly2); polys(ip+1:end)];
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



