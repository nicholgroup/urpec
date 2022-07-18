function [polys] = triangulatePoly(poly)
%triangulatePoly fractures a polygon into triangles.
%   
%[out] = triangulatePoly(poly) gives a cell array of polygons (out), where
%out{i}.x and out{i}.y are column vectors specifying the points of the new
%polygons. The input polygon (parant) is a struct with fields x
%and y, which are column vectors of the polygons.

%First, make a nice polygon.
polyin=polyshape(poly.x,poly.y);

x=polyin.Vertices(:,1);
y=polyin.Vertices(:,2);

T = triangulation(polyin);

%Standard triangulation
CL=T.ConnectivityList;
polys={};
for itri=1:size(T.ConnectivityList,1)
    xnew=x(CL(itri,:));
    ynew=y(CL(itri,:));
    polys{itri}.x=xnew;
    polys{itri}.y=ynew;
end

%Use Matlab's built in meshing functions. 

% %Create a PDE model.
% model = createpde;
% 
% %With the triangulation data as a mesh, use the geometryFromMesh function to create a geometry. Plot the geometry.
% tnodes = T.Points';
% telements = T.ConnectivityList';
% geometryFromMesh(model,tnodes,telements);
% 
% Hmin=max((max([x y])-min([x y]))/sqrt(length(x(:)))); %try to subdivide into 5 triangles.
% Hmin=max([Hmin,.1]); %minimum size of 100 nm.
% m=generateMesh(model,'Hmin',Hmin,'GeometricOrder','linear');
% %m=generateMesh(model,'GeometricOrder','linear');
% 
% 
% nx=m.Nodes(1,:);
% ny=m.Nodes(2,:);
% 
% out={};
% for im=1:size(m.Elements,2)
%     out{im}.x=nx(m.Elements([1:3],im));
%     out{im}.y=ny(m.Elements([1:3],im));
% end

end

