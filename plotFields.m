function [] = plotFields(filename)
%plotFields will plot the shapes in a fields file along with those shapes fractured into shapes with four or fewer sides.

ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0] };

if ~exist('filename','var')
    [filename, pathname]=uigetfile('*_fields.mat');
    [pathname,filename,ext] = fileparts(fullfile(pathname,filename));
    curdir=pwd;
    cd(pathname);
elseif isempty(filename)
    [filename, pathname]=uigetfile('*_fields.mat');
    [pathname,filename,ext] = fileparts(fullfile(pathname,filename));
    curdir=pwd;
    cd(pathname);
else
    curdir=pwd;
end

d=load(filename);
fields=d.fields;

if isfield(fields,'ctab')
    ctab=fields.ctab;
end

%Make a colormap based on these colors for later use.
cmap=[];
for iColor=1:length(ctab)
    cmap=[cmap; ctab{iColor}./255];
end

polygons=fields.polygons;

figure(1); clf; hold on;
colormap(cmap);

for ipoly=1:length(polygons)
    
    p=polygons(ipoly).p;
    
    x=polygons(ipoly).p(:,1);
    y=polygons(ipoly).p(:,2);
    
    [x,y,polyin]=fixPoly(x,y);
    
    nvertices=size(polyin.Vertices,1);
        
    if isempty(x)
        continue
    end
        
    %subplot(1,2,1);
    plot([x; x(1)],[y; y(1)],'color',polygons(ipoly).color./255)
    
    if ~any(nvertices==[3,4])
        %fprintf('Unsupported polygon size. Fracturing...\n');
        T = triangulation(polyin);
        
        CL=T.ConnectivityList;
        for itri=1:size(T.ConnectivityList,1)
            xnew=[x(CL(itri,:)); x(CL(itri,1))];
            ynew=[y(CL(itri,:)); y(CL(itri,1))];
            %subplot(1,2,1); hold on;
            plot(xnew,ynew,'color',polygons(ipoly).color./255);
%             subplot(1,2,2); hold on;
%             plot(xnew,ynew,'color','r');
        end
    end
    
    
end

%Add a colorbar to the plot
set(gca, 'CLim', [fields.dvalsAct(1), fields.dvalsAct(end)]);
colorbar;    
xlabel('x (microns)');
ylabel('y (microns)');
box on;
    
cd(curdir)
end

