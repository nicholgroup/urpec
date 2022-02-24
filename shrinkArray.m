function [xinds,yinds] = shrinkArray(xp,yp,p)
%shrinkArray take the subset of a spatial matrix around a polygon
%   
% shrinkArray finds the minimal subset of x values and y values in xp and
% yp (arrays of points in the x and y directions) that is needed to be
% larger than p, which is a 2-column array of x and y coordinates for a
% polygon.
% xinds and yinds are the output index values that span the polygon.
try
    indx(1)=find(xp<min(p(:,1)),1,'last');
catch
    indx(1)=1;
end

try
    indx(2)=find(xp>max(p(:,1)),1,'first');
catch
    indx(2)=length(xp);
end

try
    indy(1)=find(yp<min(p(:,2)),1,'last');
catch
    indy(1)=1;
end

try
    indy(2)=find(yp>max(p(:,2)),1,'first');
catch
    indy(2)=length(yp);
end

xinds=(indx(1):1:indx(2));
yinds=(indy(1):1:indy(2));
end

