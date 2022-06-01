function [  ] = urpec_testExpose( config )
% urpec_testExpose shows the expected developed pattern in electron beam lithography.
%
% config is an optional struct with the following optional fields:
%
%   dx: spacing in units of microns for the deconvolution. The default is
%   0.1 mcirons. It is also best to have the step size several times
%   larger than the center-center or line-line spacing in npgs. 
%   
%   targetPoints: approximate number of points for the simulation. Default
%   is 50e6.
%
%   autoRes: enables auto adjustment of the resolution for ~10min
%   computation time
%
%   file: datafile for processing. This can either be a .dxf file or a .mat
%   file. If it is a .mat file, the contets of the file should be a struct
%   array called polygons. Each element of the struct array should have
%   these fields:
%   p: an nx2 array of coordinates describing the poylgon.
%   layer: a number giving the polygon number
%   according to the convention described above.
%   Optional fields in this struct include:
%   lineType: linetype for .dc2 and NPGS
%   dose: The dose (in units to clear) of the polygon
%
%   psfFile: point-spread function file
%
%   padLen: the size with which to pad the CAD file, in units of microns.
%   The defaul is 5 microns. In general, a larger pad size will improve the
%   accuracy by accounting for long-distance proximity effects, but the
%   computation will take longer.
%
%
%
% By:
% Adina Ripin aripin@u.rochester.edu
% Elliot Connors econnors@ur.rochester.edu
% John Nichol jnich10@ur.rochester.edu
%

tic

debug=0;

if ~exist('config','var')
    config=struct();
end

config = def(config,'dx',.1);   %Grid spacing in microns. This is now affected by config.targetPoints.
config = def(config,'targetPoints',5e6);  %Target number of points for the simulation. 50e6 takes about 10 min to complete.
config = def(config,'autoRes',true);  %auto adjust the resolution
config=def(config,'file',[]); 
config=def(config,'psfFile',[]);
config=def(config,'padLen',5);

%used below. These jet-like colors are compatible with NPGS.
ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0] };

dx = config.dx;

fprintf('urpec is running...\n');

% ########## Load pattern file ##########

if isempty(config.file)
    %choose and load file
    fprintf('Select your cad file.\n')
    [filename, pathname]=uigetfile({'*.dxf';'*.mat'});
    [pathname,filename,ext] = fileparts(fullfile(pathname,filename));
else
    [pathname,filename,ext] = fileparts(config.file);

end

pathname=[pathname '\'];
filename=[filename ext];

config=def(config,'outputDir',pathname);
if config.outputDir(end)~='\'
    config.outputDir=[config.outputDir '\'];
end
    
if strmatch(ext,'.dxf')
    [lwpolylines,lwpolylayers]=dxf2coord_20(pathname,filename);
    %These are actually the layer names, not the numbers. 
    %If they do not follow the convention, then default to layer 2 and
    %fracturing.
    for i=1:length(lwpolylayers)
        if isempty(str2num(lwpolylayers{i}))
            lwpolylayers{i}='2';
        end
    end
    layerNum=str2num(cell2mat(lwpolylayers));
end

if strmatch(ext,'.mat')
    d=load([pathname filename]);
    
    if isfield(d,'fields')
        d=d.fields;
    end
    
    layerNum=[d.polygons.layer];
    lwpolylines=[ones(size(d.polygons(1).p,1),1).*1 d.polygons(1).p];
    for i=2:length(d.polygons)
        lwpolylines=[lwpolylines; [ones(size(d.polygons(i).p,1),1).*i d.polygons(i).p]];
    end
    
    if isfield(d.polygons,'dose')
        dose=[d.polygons.dose];
    else
        dose=ones(1,length(d.polygons));
    end
end

fprintf('CAD file imported.\n')

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

%defult to fracturing if the layer number extraction went bad. This can
%happen with the wrong autcad version, for example.
if length(layerNum)~=length(objects)
    layerNum=ones(1,length(objects)).*2; 
end

display(['CAD file analyzed.']);

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
   
maxXold=round(maxX,4);
minXold=round(minX,4);
maxYold=round(maxY,4);
minYold=round(minY,4);
padSize=ceil(config.padLen/dx).*dx;
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
    padSize=ceil(config.padLen/dx).*dx;
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

%make sure the number of points is odd. This is important for deconvolving
%the psf
addX=0;
if ~mod(length(xp),2)
    %maxX=maxX+dx;
    addX=1;
    %xp=[xp; xp(end)+dx]; %QHF 2021/3/30: I don't think this is correct
    xp=[xp xp(end)+dx];
end
addY=0;
if ~mod(length(yp),2)
    %maxY=maxY+dx;
    addY=1;
    yp=[yp yp(end)+dx];
end

% xp = minX:dx:maxX;
% yp = minY:dx:maxY;
[XP, YP] = meshgrid(xp, yp);

[mp, np] = size(XP);

totgridpts = length(xp)*length(yp);
polysbin = zeros(size(XP));

progressbar('Analyzing polygons');
for ar = 1:length(objects) 
    progressbar(ar/length(objects));
    p = objects{ar};
    
    subpoly = inpolygon(XP, YP, p(:,2), p(:,3));
    %polysbin = or(polysbin, subpoly);
    
    polysbin=polysbin+subpoly.*dose(ar); 
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
psfForward=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2));
psfBackscatter=1/(1+eta).*(eta/(pi*beta^2).*exp(-rpsf2./beta.^2));
%psfTot=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2)+eta/(pi*beta^2).*exp(-rpsf2./beta.^2));
%Now force the ratio of the forward and backscattered parts to be the
%correct.
psf=psfForward./sum(psfForward(:))+eta*psfBackscatter./sum(psfBackscatter(:));

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

