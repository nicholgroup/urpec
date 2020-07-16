function [  ] = urpec_v3( config )
% function [  ] = urpec_v3( config )
% Generates a proximity-effect-corrected .dxf file for electron beam 
% lithography.
% 
% The corrected file is created by deconvolving a point spread function 
% from an input .dxf or .mat pattern file.
% 
% The output file has different colors, each of which recieve a different
% dose. This function assumes that one unit in the input pattern file is one
% micron.
%
% The layer scheme is as follows. The names for all layers should be numbers.
% Layers 1 and 2 of the input file will
% both be output to layer 1 of the output file. Layer 1 will not be
% fractured, and layer 2 will be fractured. Layers 3 and 4 of the input
% file will be output to layer 2 of the output filed, etc. If the polygons
% are not fractured, the are written with an average dose. 
%
% Right now this is intended for use with NPGS.
%
% config is an optional struct with the following optional fields:
%
%   dx: spacing in units of microns for the deconvolution. The default is
%   0.01 mcirons. It is also best to have the step size several times
%   larger than the center-center or line-line spacing in npgs. 
%
%   maxIter: maximum number of iterations for the deconvolution. Default is
%   6.
%
%   dvals: doses corresponding to the layers in the output file. Default is
%   1, 1.1, 2.0 in units of the dose to clear.
%   
%   targetPoints: approximate number of points for the simulation. Default
%   is 50e6.
%
%   autoRes: enables auto adjustment of the resolution for ~10min
%   computation time
%
%   file: datafile for processing. This can either be a .dxf file or a .mat
%   file. If it is a .mat file, the contets of the file should be a struct
%   called polygons. The polygons struct should have at least these fields:
%       p: a cell array of polygons. Each element of the cell array should
%       be a nx2 array of coordinates describing the poylgon.
%       layer: an array of numbers specifying the layer of each polygon
%       according to the convention described above.
%
%   psfFile: point-spread function file
%
%
%   fracNum: maximum number of times to fracutre a shape
%
%   fracSize: minimum size for fractured shapes, in units of dx.
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
% Version history
% v2: handles different write fields and writes directly to dc2 format.
% v3: 
%   writes all doses to the same layer but with different colors. 
%   PSF improvements
%   Entirely new fracturing algorithm
%

tic

debug=0;

if ~exist('config','var')
    config=struct();
end

config = def(config,'dx',.01);   %Grid spacing in microns. This is now affected by config.targetPoints.
config = def(config,'targetPoints',50e6);  %Target number of points for the simulation. 50e6 takes about 10 min to complete.
config = def(config,'autoRes',true);  %auto adjust the resolution
config = def(config,'maxIter',6);  %max number of iterations for the deconvolution
config = def(config,'dvals',linspace(1,2.0,15));  %doses corresponding to output layers, in units of dose to clear
config=def(config,'file',[]); 
config=def(config,'psfFile',[]);
config=def(config,'fracNum',3);
config=def(config,'fracSize',4);


%used below. These jet-like colors are compatible with NPGS.
ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0] };

dx = config.dx;

fprintf('urpec is running...\n');

if isempty(config.file)
    %choose and load file
    fprintf('Select your cad file.\n')
    [filename, pathname,ext]=fileparts(uigetfile({'*.mat';'*.dxf'}));
else
    [pathname,filename,ext] = fileparts(config.file);
    pathname=[pathname '\'];
    filename=[filename ext];
end

if strmatch(ext,'.dxf')
    [lwpolylines,lwpolylayers]=dxf2coord_20(pathname,filename);
    %For later use in breaking into fields
    %These are actually the layer names, not the numbers. 
    %If they do not follow the convention, then default to layer 1 and
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
    layerNum=[d.polygons.layer];
    lwpolylines=[ones(size(d.polygons(1).p,1),1).*1 d.polygons(1).p];
    for i=2:length(d.polygons)
        lwpolylines=[lwpolylines; [ones(size(d.polygons(i).p,1),1).*i d.polygons(i).p]];
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

%make sure the number of points is odd. This is important for deconvolving
%the psf
addX=0;
if ~mod(length(xp),2)
    maxX=maxX+dx;
    addX=1;
end
addY=0;
if ~mod(length(yp),2)
    maxY=maxY+dx;
    addY=1;
end

xp = minX:dx:maxX;
yp = minY:dx:maxY;
[XP, YP] = meshgrid(xp, yp);

[mp, np] = size(XP);

totgridpts = length(xp)*length(yp);
polysbin = zeros(size(XP));


%Load the fields file. This has the locations of the field centers as well
%as the shapes of the write fields.
% if config.field
%     if isempty(config.fieldFile)
%         %choose and load file
%         fprintf('Select your fields file.\n')
%         [fname pname]=uigetfile('*.mat');
%     else
%         [pname,fname,ext] = fileparts(config.fieldFile);
%         pname=[pname '\'];
%         fname=[fname ext];
%     end
%     
%     load([pname fname]);
% else
%     fields=struct();
%     fields.center=[0 0];
% end
% 
% for i=1:length(fields)
%     fields(i).fieldShift=-fields(i).center;
% end

for ar = 1:length(objects) %EJC 5/5/2018: run time (should) scale ~linearly~ with med/sm object num
    p = objects{ar};
    
    subpoly = inpolygon(XP, YP, p(:,2), p(:,3));
    %polysbin = or(polysbin, subpoly);
    
    polysbin=polysbin+subpoly;
    
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
minSize=min([max(abs(XP(:))),max(abs(YP(:)))]);
psfRange=round(min([minSize,20]));
descr=psf.descr;

%target dx for making the psf. We will average later.
dx_t=.002;
dx_f=round(dx/dx_t);

%psf changes definition immediately below. Note: this scheme only works
%properly if the step size is larger than the 
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
    padPoints1=padPoints;
elseif xpad<0
    polysbin=padarray(polysbin,[-xpad/2,0],0,'both');
    padPoints1=padPoints-xpad/2;
end

%pad in the y direction
ypad=size(polysbin,2)-size(psf,2);
if ypad>0
    psf=padarray(psf,[0,ypad/2],0,'both');
    padPoints2=padPoints;
elseif ypad<0
    polysbin=padarray(polysbin,[0,-ypad/2],0,'both');
    padPoints2=padPoints-ypad/2; 
end

%normalize
psf=psf./sum(psf(:));

dstart=polysbin;
shape=polysbin>0; %1 inside shapes and 0 everywhere else

dose=dstart;
doseNew=shape; %initial guess at dose. Just the dose to clear everywhere
figure(555); clf; imagesc(xp,yp,polysbin);
set(gca,'YDir','norm');
title('CAD pattern');
drawnow;

%iterate on convolving the psf, and add the difference between actual dose and desired dose to the programmed dose
meanDose=0;
for i=1:config.maxIter
    %convolve with the point spread function, 
    %doseActual=ifft2(fft2(doseNew).*fft2(psf));
    doseActual=ifft2(fft2(doseNew).*fft2(psf)); 

    doseActual=real(fftshift(doseActual)); 
    %The next line is needed becase we are tryin to do FFT shift on an
    %array with an odd number of elements. 
    doseActual(2:end,2:end)=doseActual(1:end-1,1:end-1);
    doseShape=doseActual.*shape; %total only given to shapes. Excludes area outside shapes. We don't actually care about this.
    meanDose=nanmean(doseShape(:))./mean(shape(:))
    
    figure(556); clf;
    subplot(1,2,2);
    imagesc(yp,xp,doseActual);
    title(sprintf('Actual dose. Iteration %d',i));
    set(gca,'YDir','norm');
    
    doseNew=doseNew+1.2*(shape-doseShape); %Deonvolution: add the difference between the desired dose and the actual dose to doseShape, defined above
    subplot(1,2,1);
    imagesc(yp,xp,doseNew);
    title(sprintf('Programmed dose. Iteration %d',i));
    set(gca,'YDir','norm');
    
    drawnow;
    
end

dd=doseNew;
ss=shape;
%unpad arrays
doseNew=dd(padPoints1+1:end-padPoints1-1*addY,padPoints2+1:end-padPoints2-1*addX);
shape=ss(padPoints1+1:end-padPoints1-1,padPoints2+1:end-padPoints2-1);
mp=size(doseNew,1);
np=size(doseNew,2);

fprintf('done.\n')

%This is needed to not count the dose of places that get zero or NaN dose.
try
    doseNew(doseNew==0)=NaN;
    doseNew(doseNew<0)=NaN;
end

figure(557); clf;
hist(doseNew(:),config.dvals);
xlabel('Relative dose');
ylabel('Count');
drawnow;

doseStore=doseNew;

dvals=config.dvals;
nlayers=length(dvals);
dvalsl=dvals-(dvals(2)-dvals(1));
dvalsr=dvals;
layer=[];
figSize=ceil(sqrt(nlayers));
dvalsAct=[];
 

%Uncomment to save some stuff for practicing.
% d.xp=xp;
% d.yp=yp;
% d.xpold=xpold;
% d.ypold=ypold;
% d.objects=objects;
% d.doseNew=doseNew;
% d.dvals=dvals;
% %save('fractureData','d');

%Check for concave polygons. We will break them up into triangles for
%fracturing.
objectsTmp={};
layerNumTmp=[];

for ar = 1:length(objects)
    
    p=objects{ar};
    
    p=p(:,2:3);
    
    fracture=~mod(layerNum(ar),2);

    isConvex = checkConvex(p(:,1)',p(:,2)');
    if ~isConvex && fracture
        fprintf('Concave polygon to be fractured found. \n');
        x=p(:,1)';
        y=p(:,2)';
        
        %Convert polygon to polyshape type for matlab triangulation. If
        %this gives a warning, it means the concave polygon is really messed
        %up, and it's probably not worth continuing.
        warning('');
        %get rid of duplicates before doing polyshape. This is an easy way
        %to fix some bad polygons.
        A=[x' y'];
        [~, ia]=unique(A,'rows');
        ia=sort(ia);
        ia=[ia; 1];
        %A=A(ia,:);
        x=A(:,1)'; y=A(:,2)';
        polyin=polyshape(x,y);
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            figure(666); clf; plot(x,y);
            title('Bad polygon!');
            error('Bad concave polygon found. Fix it in your pattern, and try again.');
        end
        x=polyin.Vertices(:,1)';
        y=polyin.Vertices(:,2)';
        T = triangulation(polyin);
        %triplot(T);
        CL=T.ConnectivityList;
        for tt=1:length(T.ConnectivityList)
            xnew=[x(CL(tt,:)) x(CL(tt,1))];
            ynew=[y(CL(tt,:)) y(CL(tt,1))];
            znew=ones(length(xnew),1);
            objectsTmp{end+1}=[znew xnew' ynew'];
            layerNumTmp(end+1)=layerNum(ar);
        end
        
    else
        objectsTmp{end+1}=objects{ar};
        layerNumTmp(end+1)=layerNum(ar);
    end
            
end
objects=objectsTmp;
layerNum=layerNumTmp;

fprintf('Fracturing...\n');
fracCount=0;
notFracCount=0;

fracDebug=0;

subField=struct();
[XPold, YPold] = meshgrid(xpold, ypold);
fracture=0;
for ar = 1:length(objects)
    
    p=objects{ar};

    
    p=p(:,2:3);
    
    pu=unique(p,'rows');
    %the dxf importer has duplicate points.
    if length(pu)==length(p)/2-1
        p=p(1:end/2,:);
    end
    
    subpoly = inpolygon(XPold, YPold, p(:,1), p(:,2));
    shotMap=doseNew.*subpoly;
    try
        shotMap(shotMap==0)=NaN;
        shotMap(shotMap<0)=NaN;
    end
    
    fracture=~mod(layerNum(ar),2);
    
    nPolys2Frac=sum(~mod(layerNum,2));
    
    nPolysNot2Frac=sum(mod(layerNum,2));
        
    if fracture
        
        fracCount=fracCount+1;
        fprintf('Fracturing polygon %d, %d of %d \n',ar,fracCount,nPolys2Frac);
        %First figure out the variation in dose across the polygon
        maxDose=max(shotMap(:));
        minDose=min(shotMap(:));
        
        %Maximum dose variation inside of one shape. 
        dDose=dvals(2)-dvals(1);
        
        %Reduce the size of the shot map to save memory
        shape=shotMap>0;
        tt=XPold.*shape;
        maxX=max(tt(:));
        minX=min(tt(:));
        sizeX=maxX-minX;
        
        tt=YPold.*shape;
        maxY=max(tt(:));
        minY=min(tt(:));
        sizeY=maxY-minY;
        clear tt;
        
        dx=xpold(2)-xpold(1);
        dy=ypold(2)-ypold(1);
        
        %reduce the size of the shot map.
        xinds=round([(minX-xpold(1))/dx+1:(maxX-xpold(1))/dx+1]);
        yinds=round([(minY-ypold(1))/dy+1:(maxY-ypold(1))/dy+1]);
        
        %figure(555); clf; imagesc(shotMap);
        
        shotMapNew=shotMap(yinds,xinds);
        XPnew=XPold(yinds,xinds)+dx/2; xpnew=xpold(xinds)+dx/2;
        YPnew=YPold(yinds,xinds)+dy/2; ypnew=ypold(yinds)+dy/2;
        
%         figure(556); clf; imagesc(xpnew,ypnew,shotMapNew);
%         hold on; plot(p(:,1),p(:,2));
        
        minSize=config.dx*config.fracSize; %Smallest eventual polygon size
        fracNum=config.fracNum; %number of times to fracture each polygon each iteration.
        
        if (maxDose-minDose)<dDose %don't need to fracture          
            subField(ar).poly(1).x=p(:,1);
            subField(ar).poly(1).y=p(:,2);
            [mv,ind]=min(abs(dvals-squeeze(nanmean(shotMap(:)))));
            d=i; %dose color
            subField(ar).poly(1).dose=ind;
            subField(ar).poly(1).layer=ceil(layerNum(ar)/2);
        else %Need to fracture. In the future, this should be made into a function if increased complexity is desired.        
            poly=struct;
            poly(1).x=p(:,1);
            poly(1).y=p(:,2);
            poly(1).good=0;
            poly(1).layer=ceil(layerNum(ar)/2);
            [mv,ind]=min(abs(dvals-nanmean(shotMap(:))));
            poly(1).dose=ind; %dose color
            
            allGood=0;
            
            nPolys=1;
            iter=0;
            
            while ~allGood
                %Show the fracturing status. Also check for empty polygons
                %and zero-area polygons.
                if fracDebug
                    figure(777); clf; hold on;
                end
                for j=(length(poly):-1:1)
                    if isempty(poly(j).x)
                        poly(j)=[];
                    elseif polyarea(poly(j).x,poly(j).y)==0
                        poly(j)=[];
                    else
                        if fracDebug
                            plot(poly(j).x,poly(j).y,'Color',ctab{poly(j).dose}./255);
                        end
                    end
                end
                %drawnow;
                iter=iter+1;
                i=1;
                while i<=nPolys              
                    i;
                    if ~poly(i).good
                        poly(i).sizeX=max(poly(i).x)-min(poly(i).x);
                        poly(i).sizeY=max(poly(i).y)-min(poly(i).y);
                    
%                       %In preparation for fracturing, reduce the size of the shot map.
                        maxX=max(poly(i).x);
                        maxY=max(poly(i).y);
                        minX=min(poly(i).x);
                        minY=min(poly(i).y);                        
                        xinds=round([(minX-xpold(1))/dx+1:(maxX-xpold(1))/dx]);
                        yinds=round([(minY-ypold(1))/dy+1:(maxY-ypold(1))/dy]);                                                
                        shotMapNew=shotMap(yinds,xinds);
                        XPnew=XPold(yinds,xinds)+dx/2; xpnew=xpold(xinds)+dx/2;
                        YPnew=YPold(yinds,xinds)+dy/2; ypnew=ypold(yinds)+dy/2;
                        
                        xline=squeeze(nanmean(shotMapNew,1));
                        yline=squeeze(nanmean(shotMapNew,2));
                        xdiff=max(xline)-min(xline);
                        ydiff=max(yline)-min(yline);
                        
                        shouldFracX=(xdiff>dDose/(2+iter*.5));
                        if isempty(shouldFracX)
                            shouldFracX=0;
                        end
                        shouldFracY=(ydiff>dDose/(2+iter*.5));
                        if isempty(shouldFracY)
                            shouldFracY=0;
                        end
                        canFracX=(poly(i).sizeX/fracNum>minSize);
                        canFracY=(poly(i).sizeY/fracNum>minSize);
                        
%                         if isempty(canFracX) || isempty(shouldFracX) || isempty(canFracY) ||isempty(shouldFracY)
%                             error('Fracturing error. Try decreasing the step size')
%                         end
                        
                        %Actually do the fracturing
                        if canFracX || canFracY
                            polys2Add=DIVIDEXY(poly(i),canFracX*shouldFracX*(fracNum-1)+1,canFracY*shouldFracY*(fracNum-1)+1);
                            polys2Add=polys2Add(:);
                        else
                            polys2Add={};
                        end
                        
                        %Check for empty polys and get rid of them
                        %Also check for polys of zero area and get rid of
                        %them.
                        for j=(length(polys2Add):-1:1)
                            
                            if isempty(polys2Add{j}) || any(size(polys2Add{j}.x)~=size(polys2Add{j}.y))
                                polys2Add(j)=[];
                            end
                            
                            try 
                                if polyarea(polys2Add{j}.x,polys2Add{j}.y)<1e-5
                                    polys2Add(j)=[];
                                end
                            end
                            
                        end

                        %If after all this we didn't actually fracture anything,
                        %then skip and move on.
                        if length(polys2Add)==1
                            polys2Add={};
                            i=i+1;
                        end
                        
                        length(polys2Add);
                        
                        %Go through the new polygons and find out if
                        %they're good.
                        for j=1:length(polys2Add)
                            %plot(polys2Add{j}.x,polys2Add{j}.y)
                            %Close the polys.
                           
                            if polys2Add{j}.x(1)~=polys2Add{j}.x(end) || polys2Add{j}.y(1)~=polys2Add{j}.y(end)
                                polys2Add{j}.x= [polys2Add{j}.x polys2Add{j}.x(1)];
                                polys2Add{j}.y= [polys2Add{j}.y polys2Add{j}.y(1)];
                            end
                            polys2Add{j}.x= polys2Add{j}.x';
                            polys2Add{j}.y= polys2Add{j}.y';                            
                            polys2Add{j}.good=0;
                            polys2Add{j}.sizeX=max(poly(i).x)-min(poly(i).x);
                            polys2Add{j}.sizeY=max(poly(i).y)-min(poly(i).y);
                            polys2Add{j}.layer=ceil(layerNum(ar)/2);
                            
                            subpoly = inpolygon(XPnew, YPnew, polys2Add{j}.x(:), polys2Add{j}.y(:));
                            sm=shotMapNew.*subpoly;
                            sm(sm==0)=NaN;
                            maxDose=max(sm(:));
                            minDose=min(sm(:));
                            
                            [mv,ind]=min(abs(dvals-nanmean(sm(:))));
                            polys2Add{j}.dose= ind;

                            if (maxDose-minDose)<dDose %don't need to fracture
                                polys2Add{j}.good=1;                                
                            end
                        end
                        
                        %Add the new polygons back to the total list of polygons
                        if ~isempty(polys2Add)
                            if iter==1 
                                for j=1:length(polys2Add)
                                    poly(j)=polys2Add{j};
                                end
                            else
                                polysEnd=poly(i+1:end);
                                %Add the new polygons in the right place.
                                poly=[poly(1:i-1) polys2Add{:} polysEnd];
                            end
                        end
                        
                        %Update the tot polygon number and counter (i)
                        nPolys=length(poly);
                        i=i+length(polys2Add);
                        
                    else %if we didn't need to fracture it, go to the next polygon.
                        i=i+1;
                    end
                    
                    %If we can't fracture anymore, call it good.
                    if (~canFracX && ~canFracY) 
                        poly(i).good=1;
                    end
                end
                
                if iter==fracNum
                    for j=1:length(poly)
                        poly(j).good=1;
                    end
                end
                allGood=all([poly.good]);    
                
            end
            
            subField(ar).poly=poly;                      
        end
        
        dvalsAct=dvals;
        
    else %no fracturing
        
        notFracCount=notFracCount+1;
        fprintf('Averaging polygon %d, %d of %d \n',ar,notFracCount,nPolysNot2Frac);
        
        subField(ar).poly(1).x=p(:,1);
        subField(ar).poly(1).y=p(:,2);
        [mv,ind]=min(abs(dvals-squeeze(nanmean(shotMap(:)))));
        subField(ar).poly(1).dose=ind;
        subField(ar).poly(1).layer=ceil(layerNum(ar)/2);
        
        dvalsAct=dvals;
        
    end

end

fprintf('Fracturing complete. \n')

%save the final files
fields=struct();

j=1; %residual counter left over from when we had a fields struct.
outputFileName=[pathname filename(1:end-4) '_' descr '_' num2str(j) '.dxf'];
dc2FileName=[pathname filename(1:end-4) '_' descr '_' num2str(j) '.dc2'];

fprintf('Exporting to %s...\n',outputFileName);

fields(j).cadFile=[filename(1:end-4) '_' descr '_' num2str(j) '.dc2'];

try
    FID = dxf_open(outputFileName);
catch
    fprintf('Close the dxf file, and then type dbcont and press enter. \n')
    keyboard;
    FID = dxf_open(outputFileName);
end
        
polygons=struct();
polygons(1)=[];

%Write in order of objects.
for ar=1:length(subField)
    
    for b=1:length(subField(ar).poly)
        if ~isempty(subField(ar).poly(b).x) %Needed because sometimes our fracturing algorithm generates empty polygons.
            i=subField(ar).poly(b).dose;
            FID=dxf_set(FID,'Color',ctab{i}./255,'Layer',subField(ar).poly(b).layer); 
            X=subField(ar).poly(b).x;
            Y=subField(ar).poly(b).y;
            Z=X.*0;
            dxf_polyline(FID,X,Y,Z);
            %plot(X,Y);
            
            polygons(end+1).p=[X Y];
            polygons(end).color=ctab{i}; %JMN 2019/10/12
            polygons(end).layer=subField(ar).poly(b).layer;
            polygons(end).lineType=1;
        end
        
    end
    
end


dxf_close(FID);

fprintf('Exporting to %s...\n',dc2FileName);

dc2write(polygons,dc2FileName);

%Save doses here
doseFileName=[pathname filename(1:end-4) '_' descr '_' num2str(j) '.txt'];

fileID = fopen(doseFileName,'w');
fprintf(fileID,'%3.3f \r\n',dvalsAct);
fclose(fileID);

%Save all of the information in a .mat file for later.
fields(j).doseFile=[filename(1:end-4) '_' descr '_' num2str(j) '.txt'];
fields(j).dvalsAct=dvalsAct;
fields(j).polygons=polygons;
fprintf('Finished exporting.\n');

fieldsFileName=[pathname filename(1:end-4) '_' descr '_fields.mat'];
save(fieldsFileName,'fields');

fprintf('urpec is finished.\n')

toc

end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
end

