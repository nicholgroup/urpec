function [] = dc2write(polygons,fname)
%[] = dc2write(polygons,fname)
%   Write polygons into a dc2 file
%
%   fname: filename
%   polygons: struct array with the following fields
%       p: polygon
%       layer: layer number
%       color: color in rgb format
%       line: linetype
%
% the combination \r\n is somehow essential for NPGS.


pblock=[];
minX=0;
minY=0;
maxX=0;
maxY=0;
layers=[];

for j=1:length(polygons)
    layers(j)=polygons(j).layer;
end

uniqueLayers=unique(layers);
if uniqueLayers(1)~=0
    uniqueLayers=[0 uniqueLayers];
end

for k=1:length(uniqueLayers)
    inds=(layers==uniqueLayers(k));
    if ~(any(inds))
        continue
    end
    %layer line
    pblock=[pblock sprintf('21 %d 0 0 0 0\r\n',k)];
    %add all polygons
    for j=1:length(polygons)
        if ~inds(j)
            continue
        end
        npoints=size(polygons(j).p,1);

        %pblock=[pblock sprintf('1 %d 8 0 %d 10 0 1 0 %d %d %d \r\n',npoints,polygons(j).lineType,polygons(j).color(1),polygons(j).color(2),polygons(j).color(3))];
        
        %New element header compatible with NPGS
        pblock=[pblock sprintf('1 %d 8 0 %d 10 0 1 0 %d %d %d 0 1 \r\n',npoints,polygons(j).lineType,polygons(j).color(1),polygons(j).color(2),polygons(j).color(3))];

        for i=1:size(polygons(j).p,1)
            x=polygons(j).p(i,1)*8;
            y=polygons(j).p(i,2)*8;
            pblock=[pblock sprintf('%.4f %.4f 0\r\n',x,y)];
            if x<minX
                minX=x;
            end
            if y<minY
                minY=y;
            end
            if x>maxX
                maxX=x;
            end
            if y>maxY
                maxY=y;
            end
        end
    end
end

lx=maxX-minX;
ly=maxY-minY;

hblock=[];
% hblock=sprintf('%.4f %.4f %.4f %.4f 0 0\r\n',minX,minY,lx,ly);
% hblock=[hblock sprintf('20 0 0 0 0 0\r\n8.000000, 0.800000\r\n8.000000\r\n*\r\n23 %d 0 0 0 0\r\n\r\n',length(uniqueLayers)+1)];

%New header comaptible with NPGS
hblock=sprintf('%.4f %.4f %.4f %.4f 0 -0.0000 0.0000\r\n',minX,minY,lx,ly);
hblock=[hblock sprintf(['42 20 0 0 0 0\r\n8.000000\r\n8.000000, 0.800000\r\n3.000000\r\n16\r\n0.000000\r\n0.000000\r\n1.000000\r\n1.000000\r\n'...
    '0\r\n0\r\n1\r\nSIMPLEX2.VFN\r\n'...
    '0.00000000 0.00000000 0.00000000 0.00000000\r\n'...
    '0.00000000 0.00000000 0.00000000 0.00000000\r\n'...
    '0.00000000 0.00000000 0.00000000 0.00000000\r\n'...
    '0.00000000 0.00000000 0.00000000 0.00000000\r\n'...
    '1 0 0 0 0 0 0 0\r\n'...
    '0.00000000 0.00000000 0.00000000 0.00000000\r\n'...
    '0.00000000 0.00000000 0.00000000 0.00000000\r\n'...
    '; DesignCAD Drawing Comments /w '';'' as 1st char.\r\n'...
    '23 %d 0 0 0 0\r\n'...
    'Do Not Use\r\n'],length(uniqueLayers)+1)];


for k=1:length(uniqueLayers)
    hblock=[hblock sprintf('%d\r\n',uniqueLayers(k))];
end

writeBlock=[hblock pblock];

fid=fopen(fname,'w');

fprintf(fid, writeBlock);

fclose(fid);

end

