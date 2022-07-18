function [polys] = fixPoly2(poly)
%fixPoly2 tries to fix bad polygons by removing duplicate vertices and
%simplifying the shape if possible. fixPoly2 will respect polygons with
%multiple regions.
%
% poly is a struct with the fields x and y. The output argument is a struct
% array of polygons.

x=poly.x(:);
y=poly.y(:);

polys=poly;

%Clean up the shape.
%Get rid of duplicates before doing polyshape. This is an easy way
%to fix some bad polygons.
A=[x y];
[C, ia,ic]=unique(A,'rows');
ia=sort(ia);
ia=[ia; ];
A=A(ia,:);
xx=A(:,1); yy=A(:,2);

%try to make a nice polygon out of the coordinates. This will try
%to simplify any strange shapes and remove duplicates.
polyin=polyshape(xx,yy);
polysTmp=regions(polyin);
[warnMsg, warnId] = lastwarn;

if length(polysTmp)==0
    polys={};
end

%If we got any warnings, do it again to see if they persist
for ip=1:length(polysTmp)
    if ~isempty(warnMsg)
        lastwarn('');

        polyin=polyshape(polysTmp(ip).Vertices(:,1)',polysTmp(ip).Vertices(:,2)');
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            figure(666); clf; plot(poly.x,poly.y);
            title('Bad polygon!');
            error('Bad polygon found. Fix it in your pattern.');
        else
            x=polyin.Vertices(:,1);
            y=polyin.Vertices(:,2);            
            polys(ip).x=x;
            polys(ip).y=y;
            polys(ip).p=[x y];
        end
    else
        x=polysTmp(ip).Vertices(:,1);
        y=polysTmp(ip).Vertices(:,2);
        polys(ip).x=x;
        polys(ip).y=y;
        polys(ip).p=[x y];

    end
end

end

