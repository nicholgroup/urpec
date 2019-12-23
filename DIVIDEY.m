function [polygon1 polygon2]=DIVIDEY(polygon,Y)
%version 1.4 (2.44 KB) by Ayad Al-Rumaithi
%https://www.mathworks.com/matlabcentral/fileexchange/71635-divide-polygon
polygon1=[];
polygon2=[];
if not(isempty(polygon))
    m=length(polygon.y);
    c1=0;
    c2=0;
    for i=1:1:m
        j=i+1;
        if i==m
            j=1;
        end
        if polygon.y(i)<=Y
            c1=c1+1;
            polygon1.x(c1)=polygon.x(i);
            polygon1.y(c1)=polygon.y(i);
        end
        if polygon.y(i)>=Y
            c2=c2+1;
            polygon2.x(c2)=polygon.x(i);
            polygon2.y(c2)=polygon.y(i);
        end
        if (polygon.y(i)>Y && polygon.y(j)<Y) || (polygon.y(i)<Y && polygon.y(j)>Y)
            c1=c1+1;
            polygon1.x(c1)=polygon.x(j)+(polygon.x(i)-polygon.x(j))/(polygon.y(i)-polygon.y(j))*(Y-polygon.y(j));
            polygon1.y(c1)=Y;
            c2=c2+1;
            polygon2.x(c2)=polygon.x(j)+(polygon.x(i)-polygon.x(j))/(polygon.y(i)-polygon.y(j))*(Y-polygon.y(j));
            polygon2.y(c2)=Y;
        end
        
        
        
    end
end
