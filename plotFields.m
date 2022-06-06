function [] = plotFields(filename)
%plotFields will plot the shapes in a fields file along with those shapes fractured into shapes with four or fewer sides.

warning('off');

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

fprintf('Loading from %s. \n',filename);

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

plots=struct();
plots.X=[];
plots.Y=[];
plots(2:length(ctab))=plots(1);

progressbar('Analyzing polygons');
for ipoly=1:length(polygons)
    progressbar(ipoly/length(polygons));
    
    p=polygons(ipoly).p;
    dose=polygons(ipoly).dose;
    
    x=polygons(ipoly).p(:,1);
    y=polygons(ipoly).p(:,2);
    
    [x,y,polyin]=fixPoly(x,y);
    
    nvertices=size(polyin.Vertices,1);
        
    if isempty(x)
        continue
    end
    
    if ~any(nvertices==[3,4]) %Need to fracture the shape
        
        parent=struct;
        parent.x=x;
        parent.y=y;
        
        polys=divideAndTriangulatePoly(parent);
        
        for ip=1:length(polys)
            xnew=[polys{ip}.x(:); polys{ip}.x(1)];
            ynew=[polys{ip}.y(:); polys{ip}.y(1)];
            
            %plot(xnew,ynew,'color',polygons(ipoly).color./255);
            %Manipulating the plot is faster if we combine as much as
            %possible into one line.
            plots(dose).X=[plots(dose).X; xnew(:); NaN];
            plots(dose).Y=[plots(dose).Y; ynew(:); NaN];
        end


    else
        %plot([x; x(1)],[y; y(1)],'color',polygons(ipoly).color./255);
        %Manipulating the plot is faster if we combine as much as
        %possible into one line.
        plots(dose).X=[plots(dose).X; x; x(1); NaN];
        plots(dose).Y=[plots(dose).Y; y; y(1); NaN];

    end
    
end

%Plot each color as one line.
for ip=1:length(plots)
    plot(plots(ip).X,plots(ip).Y,'color',ctab{ip}./255)
end

%Add a colorbar to the plot
set(gca, 'CLim', [fields.dvalsAct(1), fields.dvalsAct(end)]);
colorbar;    
xlabel('x (microns)');
ylabel('y (microns)');
box on;
    
cd(curdir);

end

