function PXY=DIVIDEXY(polygon,NX,NY)
%version 1.4 (2.44 KB) by Ayad Al-Rumaithi
%https://www.mathworks.com/matlabcentral/fileexchange/71635-divide-polygon
%Input
%polygon: a structure consist of polygon.x (vector of x-coordinates) and polygon.y (vector of y-coordinates)
%NX: Number of divisions in x direction
%NY: Number of divisions in y direction
%Output
%PXY: a cell array where PX{i,j}.x and PX{i,j} are vectors of x and y coordinates of new polygon in (i,j) grid position
DX=(max(polygon.x)-min(polygon.x))/NX;
DY=(max(polygon.y)-min(polygon.y))/NY;
i=0;
P=polygon;
for X=min(polygon.x)+DX:DX:max(polygon.x)-DX
    i=i+1;
    [PX{i}, P]=DIVIDEX(P,X);
    
end
PX{NX}=P;
for i=1:1:NX
    j=0;
    for Y=min(polygon.y)+DY:DY:max(polygon.y)-DY
        j=j+1;
        [PXY{i,j}, PX{i}]=DIVIDEY(PX{i},Y);
    end
    PXY{i,NY}=PX{i};
end
end
