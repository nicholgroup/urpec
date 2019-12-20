function [polygon1 polygon2]=DIVIDEX(polygon,X)
%version 1.4 (2.44 KB) by Ayad Al-Rumaithi
%https://www.mathworks.com/matlabcentral/fileexchange/71635-divide-polygon
polygon1=[];
polygon2=[];
if not(isempty(polygon))
    m=length(polygon.x);
    c1=0;
    c2=0;
    for i=1:1:m
        j=i+1;
        if i==m
            j=1;
        end
        if polygon.x(i)<=X
            c1=c1+1;
            polygon1.x(c1)=polygon.x(i);
            polygon1.y(c1)=polygon.y(i);
        end
        if polygon.x(i)>=X
            c2=c2+1;
            polygon2.x(c2)=polygon.x(i);
            polygon2.y(c2)=polygon.y(i);
        end
        if (polygon.x(i)>X && polygon.x(j)<X) || (polygon.x(i)<X && polygon.x(j)>X)
            c1=c1+1;
            polygon1.x(c1)=X;
            polygon1.y(c1)=polygon.y(j)+(polygon.y(i)-polygon.y(j))/(polygon.x(i)-polygon.x(j))*(X-polygon.x(j));
            c2=c2+1;
            polygon2.x(c2)=X;
            polygon2.y(c2)=polygon.y(j)+(polygon.y(i)-polygon.y(j))/(polygon.x(i)-polygon.x(j))*(X-polygon.x(j));
        end
        
        
        
    end
