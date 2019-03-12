function P = simplify_polygon(P)
% function P = simplify_polygon(P)
% Removes unnecessary vertices from a polygon
% Shameless stolen from the website below.
% https://blogs.mathworks.com/steve/2012/08/28/wrapping-up-the-analysis-of-cody-solutions/
  try
  diag(sum(abs(diff(P)),2)) \ diff(P);
  P(any(ans - circshift(ans,1),2),:);
  P = vertcat(ans, ans(1,:));
  end
  
end


