function [x,y,polyin] = fixPoly(x,y)
%fixPoly tries to fix bad polygons by removing duplicate vertices and
%simplifying the shape if possible.
%
% x and y are row vectors of the polygon vertices

flipx=0; flipy=0;

if size(x,1)>1
    x=x';
    flipx=1;
end

if size(y,1)>1
    y=y';
    flipy=1;
end

warning('');

%Clean up the shape.
%Get rid of duplicates before doing polyshape. This is an easy way
%to fix some bad polygons.
A=[x' y'];
[~, ia]=unique(A,'rows');
ia=sort(ia);
ia=[ia; 1];
%A=A(ia,:);
x=A(:,1)'; y=A(:,2)';

%try to make a nice polygon out of the coordinates. This will try
%to simplify any strange shapes.
polyin=polyshape(x,y);
[warnMsg, warnId] = lastwarn;
lastwarn('');
%If polyshape gives a warning, it simpliied the shape. 
% Do it again to see if it fixed it the first time.
polyin=polyshape(polyin.Vertices(:,1)',polyin.Vertices(:,2)');
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    figure(666); clf; plot(x,y);
    title('Bad polygon!');
    error('Bad polygon found. Fix it in your pattern.');
end
if flipx
    x=polyin.Vertices(:,1);
else
    x=polyin.Vertices(:,1)'; 
end

if flipy
    y=polyin.Vertices(:,2);
else
    y=polyin.Vertices(:,2)';   
end


end

