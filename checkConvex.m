function isConvex = checkConvex(px,py)
% Given a set of points determine if they form a convex polygon
% Inputs 
%  px, py: x and y coordinates of the vertices for a given polygon

% Output 
%  isConvex: 1 or 0 for polygon being convex or concave
% https://matlabgeeks.com/tips-tutorials/computational-geometry/check-convexity-of-polygon/

numPoints = numel(px);
if numPoints < 4
    isConvex = true;
    return
end

% check to see if polygon is closed
if abs(px(end)-px(1)) < eps && abs(py(end)-py(1)) < eps
    px(end) = [];
    py(end) = [];
end
 numPoints = numel(px);

% 
% can determine if the polygon is convex based on the direction the angles
% turn.  If all angles are the same direction, it is convex.
v1 = [px(1) - px(end), py(1) - py(end)];
v2 = [px(2) - px(1), py(2) - py(1)];
signPoly = sign(det([v1; v2]));

% check subsequent vertices
for k = 2:numPoints-1
    v1 = v2;
    v2 = [px(k+1) - px(k), py(k+1) - py(k)]; 
    curr_signPoly = sign(det([v1; v2]));
    % check that the signs match
    if not (isequal(curr_signPoly, signPoly))
        isConvex = false;
        return
    end
end
% check the last vectors
v1 = v2;
v2 = [px(1) - px(end), py(1) - py(end)];
curr_signPoly = sign(det([v1; v2]));
if not (isequal(curr_signPoly, signPoly))
    isConvex = false;
else
    isConvex = true;
end