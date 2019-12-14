function [  ] = urpec_testExpose( config )
% function [  ] = urpec_testExpose( config )
% Shows the expected developed pattern in electron beam lithography.
%
% Right now this is intended for use with NPGS.
%
% config is an optional struct with the following optional fields:
%
%   dx: spacing in units of microns for the deconvolution. The default is
%   0.01 mcirons
%
%   subfieldSize: size of subfields in fracturing. Default is 50 points. This
%   is necessary for NPGS because the maximum polygon size in designCAD is
%   200 points.
%
%   maxIter: maximum number of iterations for the deconvolution. Default is
%   6.
%
%   dvals: doses corresponding to the layers in the output file. Default is
%   1, 1.1, 1.9 in units of the dose to clear.
%   
%   windowVal: smoothing distance for the dose modulation. Units are
%   approximately the grid spacing. Default is 10.
%
%   targetPoints: approximate number of points for the simulation. Default
%   is 50e6.
%
%   autoRes: enables auto adjustment of the resolution for ~10min
%   computation time
%
%   File: datafile for processing
%
%   psfFile: point-spread function file
%
%   layer: can be true or false. Choose true to respect layers in the CAD file, or
%   false not to. This is most useful if you need to break your write into
%   different fields. Default is false.
%
%   layerShift: a matrix describing an optional coordinate shift for each
%   layer. This is useful to shift the different field centers above to
%   (0,0). Can only be used in layer is 'true.' The ordering of this matrix
%   must also have the same ordering as the layers. This will break if any
%   layers are absent.
%
%   doses: an array describing the doses (relative to the DCT) for each
%   layer. Ther should be one element per layer.
%
% call this function without any arguments, or via
% urpec(struct('dx',0.005, 'subfieldSize',20,'maxIter',6,'dvals',[1:.2:2.4]))
% for example
%
%
% By:
% Adina Ripin aripin@u.rochester.edu
% Elliot Connors econnors@ur.rochester.edu
% John Nichol jnich10@ur.rochester.edu
%

tic

addpath('dxflib');

debug=0;

if ~exist('config','var')
    config=struct();
end

config = def(config,'dx',.01);   %Grid spacing in microns. This is now affected by config.targetPoints.
config = def(config,'targetPoints',50e6);  %Target number of points for the simulation. 50e6 takes about 10 min to complete.
config = def(config,'autoRes',true);  %auto adjust the resolution
config = def(config,'subfieldSize',50);  %subfield size in microns
config = def(config,'maxIter',6);  %max number of iterations for the deconvolution
config = def(config,'dvals',[1:.1:2.4]);  %doses corresponding to output layers, in units of dose to clear
config = def(config,'windowVal',10);  %Smoothing factor for the dose modultion. Default is 10. Units are approximately the grid spacing.
config=def(config,'file',[]); 
config=def(config,'psfFile',[]);
config=def(config,'layer',false);
config=def(config,'layerShift',[]);


dx = config.dx;
subfieldSize=config.subfieldSize;

fprintf('urpec is running...\n');

if isempty(config.file)
    %choose and load file
    fprintf('Select your dxf file.\n')
    [filename pathname]=uigetfile('*.dxf');
else
    [pathname,filename,ext] = fileparts(config.file);
    pathname=[pathname '\'];
    filename=[filename ext];
end

[lwpolylines,lwpolylayers]=dxf2coord_20(pathname,filename);
fields=struct();
layerNum=str2num(cell2mat(lwpolylayers));

fprintf('dxf CAD file imported.\n')

%splitting each set of points into its own object
object_num=max(lwpolylines(:,1)); % number of polygons in CAD file
shapes = lwpolylines;
objects = cell(1,object_num); % cell array of polygons in CAD file
len = size(shapes,1); % total number of polylines

% populating 'objects' (cell array of polygons in CAD file)
count = 1;
obj_num_ind_start = 1;
for obj = 1:len % loop over each polyline
    c = lwpolylines(obj,1); % c = object number
    if c ~= count
        objects{count} = lwpolylines(obj_num_ind_start:(obj-1), 1:3); %ID, x, y
        obj_num_ind_start = obj;
        count = count + 1;
        if count == object_num
            objects{count} = lwpolylines(obj_num_ind_start:(len), 1:3);
        end
    end
    
end
if object_num <= 1
    objects{count} = lwpolylines(:, 1:3);
end

%grouping objects of the same size together
areas = cell(1, object_num);
for ar = 1:object_num
    curr_obj = objects{ar};
    areas{ar} = polyarea(curr_obj(:,2), curr_obj(:,3));
end

display(['dxf CAD file analyzed.']);

if length(objects)>1
    medall=vertcat(objects{1},objects{2});
else
    medall = objects{1};
end
for i = 3:object_num
    medall = vertcat(medall, objects{i});
end

maxX = max(medall(:,2));
maxY = max(medall(:,3));
minX = min(medall(:,2));
minY = min(medall(:,3));

x = minX:dx:maxX;
y = minY:dx:maxY;
[X,Y] = meshgrid(x, y);
[m,n] = size(X);


fprintf(['Creating 2D binary grid spanning medium feature write field (spacing = ', num2str(dx), '). This can take a few minutes...\n']);
%polygon distribution - creating binary grid of x,y points that are 1 if
%in a polygon or 0 if not

%Make the simulation area bigger by 5 microns to account for proximity effects
maxXold=maxX;
minXold=minX;
maxYold=maxY;
minYold=minY;
padSize=ceil(5/dx).*dx;
padPoints=padSize/dx;
maxX=maxXold+padSize;
minX=minXold-padSize;
maxY=maxYold+padSize;
minY=minYold-padSize;

xpold = minXold:dx:maxXold;
ypold = minYold:dx:maxYold;

xp = minX:dx:maxX;
yp = minY:dx:maxY;

totPoints=length(xp)*length(yp);
fprintf('There are %2.0e points. \n',totPoints);
%Make sure the grid size is appropriate
if config.autoRes && (totPoints<.8*config.targetPoints || totPoints>1.2*config.targetPoints)
    expand=ceil(log2(sqrt(totPoints/config.targetPoints)));
    dx=dx*2^(expand);
    fprintf('Resetting the resolution to %3.4f.\n',dx);
    padSize=ceil(5/dx).*dx;
    padPoints=padSize/dx;
    maxX=maxXold+padSize;
    minX=minXold-padSize;
    maxY=maxYold+padSize;
    minY=minYold-padSize;
    xp = minX:dx:maxX;
    yp = minY:dx:maxY;
    xpold = minXold:dx:maxXold;
    ypold = minYold:dx:maxYold;
    
    fprintf('There are now %2.0e points.\n',length(xp)*length(yp));
end

%make sure the number of points is odd. This is important for deconvolving the psf
if ~mod(length(xp),2)
    maxX=maxX+dx;
end

if ~mod(length(yp),2)
    maxY=maxY+dx;
end

xp = minX:dx:maxX;
yp = minY:dx:maxY;
[XP, YP] = meshgrid(xp, yp);

[mp, np] = size(XP);

totgridpts = length(xp)*length(yp);
polysbin = zeros(size(XP));

if ~config.layer
    layerNum=layerNum.*0+1;
end

layerNums=unique(layerNum);
for i=1:length(layerNums)
    fields(i).polysbin=zeros(size(XP));
    fields(i).polys=0;
    fields(i).layerNum=layerNums(i);
    try
        fields(i).layerShift=config.layerShift(layerNums(i),:);
    catch
        fields(i).layerShift=[0,0];
    end
    
    try
        fields(i).dose=config.doses(i);
    catch
        fields(i).dose=1;       
    end
        
end
fprintf('Found %d layers. \n',length(layerNums));


for ar = 1:length(objects) %EJC 5/5/2018: run time (should) scale ~linearly~ with med/sm object num
    p = objects{ar};
    subpoly = inpolygon(XP, YP, p(:,2), p(:,3));
    %polysbin = or(polysbin, subpoly);
    
    polysbin=polysbin+subpoly*fields(layerNum(ar)).dose;
    fields(layerNum(ar)).polysbin= fields(layerNum(ar)).polysbin+subpoly;
    fields(layerNum(ar)).polys=fields(layerNum(ar)).polys+1;
end

[xpts ypts] = size(polysbin);

fprintf('done analyzing file.\n')

if isempty(config.psfFile)
    fprintf('Select point spread function file.\n')
    load(uigetfile('*PSF*'));
else
    load(config.psfFile);
end

fprintf('Deconvolving psf...');

eta=psf.eta;
alpha=psf.alpha;
beta=psf.beta;
psfRange=psf.range;
descr=psf.descr;

%psf changes definition immediately below
npsf=round(psfRange./dx);
psf=zeros(round(psfRange./dx));
[xpsf ypsf]=meshgrid([-npsf:1:npsf],[-npsf:1:npsf]);
xpsf=xpsf.*dx;
ypsf=ypsf.*dx;
rpsf2=xpsf.^2+ypsf.^2;
psf=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2)+eta/(pi*beta^2).*exp(-rpsf2./beta.^2));

%Zero pad to at least 10um x 10 um;
%pad in the x direction
xpad=size(polysbin,1)-size(psf,1);
if xpad>0
    psf=padarray(psf,[xpad/2,0],0,'both');
elseif xpad<0
    polysbin=padarray(polysbin,[-xpad/2,0],0,'both');
    for i=1:length(layerNums)
        fields(i).polysbin=padarray(fields(i).polysbin,[-xpad/2,0],0,'both');
    end
end

%pad in the y direction
ypad=size(polysbin,2)-size(psf,2);
if ypad>0
    psf=padarray(psf,[0,ypad/2],0,'both');    
elseif ypad<0
    polysbin=padarray(polysbin,[0,-ypad/2],0,'both');
    for i=1:length(layerNums)
        fields(i).polysbin=padarray(fields(i).polysbin,[-xpad/2,0],0,'both');
    end
end

%normalize
psf=psf./sum(psf(:));

doseNew=polysbin; %initial guess at dose. Just the dose to clear everywhere
figure(555); clf; imagesc(xp,yp,polysbin);
set(gca,'YDir','norm');
title('CAD pattern');
drawnow;

%iterate on convolving the psf, and add the difference between actual dose and desired dose to the programmed dose
doseActual=ifft2(fft2(doseNew).*fft2(psf)); %use window to avoid ringing
doseActual=real(fftshift(doseActual));
 
%These are the parts that develop.
pattern=doseActual>1;

figure(556); clf;
subplot(1,2,1);
imagesc(yp,xp,polysbin);
title('CAD pattern');
set(gca,'YDir','norm');


subplot(1,2,2);
imagesc(yp,xp,pattern);
title('Developed pattern');
set(gca,'YDir','norm');

fprintf('urpec is finished.\n')

toc

end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
end

