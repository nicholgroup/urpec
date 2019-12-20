function [  ] = urpec_v2( config )
% function [  ] = urpec( config )
% Generates a proximity-effect-corrected .dxf file for electron beam lithography.
% The corrected file is created by deconvolving a point spread function from an input .dxf pattern file.
% The output file has different layers, each of which should receive a different dose.
% This function assumes that one unit in the input pattern file is one
% micron.
%
% Patterns are output in the same layer in which they were input. The dose
% is modulated according to the color.
%
% Right now this is intended for use with NPGS.
%
% config is an optional struct with the following optional fields:
%
%   dx: spacing in units of microns for the deconvolution. The default is
%   0.01 mcirons. It is also best to have the step size several times
%   larger than the center-center or line-line spacing in npgs. 
%
%   fracture: turns fracturing on or off. Not immplemented yet.
%   
%
%   maxIter: maximum number of iterations for the deconvolution. Default is
%   6.
%
%   dvals: doses corresponding to the layers in the output file. Default is
%   1, 1.1, 2.5 in units of the dose to clear.
%   
%   targetPoints: approximate number of points for the simulation. Default
%   is 50e6.
%
%   autoRes: enables auto adjustment of the resolution for ~10min
%   computation time
%
%   file: datafile for processing
%
%   psfFile: point-spread function file
%
%   field: can be true or false. Choose true to if the write should be
%   broken into different write fields.
%
%   fieldFile: path to a file describing the write fields. The file should
%   have a struct array called fields, with the following fields
%       box: Bounding polygon for the write field
%      center: Stage coordinates for the field center
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
% v3: writes all doses to the same layer but with different colors. This is
% useful for controlling the write order. The writing order is now the same
% as the polygon save order.  Removed windowVal.
%

tic

debug=0;

if ~exist('config','var')
    config=struct();
end

config = def(config,'dx',.01);   %Grid spacing in microns. This is now affected by config.targetPoints.
config = def(config,'targetPoints',50e6);  %Target number of points for the simulation. 50e6 takes about 10 min to complete.
config = def(config,'autoRes',true);  %auto adjust the resolution
config = def(config,'subfieldSize',50);  %subfield size in microns
config = def(config,'maxIter',6);  %max number of iterations for the deconvolution
config = def(config,'dvals',linspace(1,2.4,15));  %doses corresponding to output layers, in units of dose to clear
config=def(config,'file',[]); 
config=def(config,'psfFile',[]);
config=def(config,'field',false);
config=def(config,'fieldFile',[]);




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
%For later use in breaking into fields
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
if config.field
    if isempty(config.fieldFile)
        %choose and load file
        fprintf('Select your fields file.\n')
        [fname pname]=uigetfile('*.mat');
    else
        [pname,fname,ext] = fileparts(config.fieldFile);
        pname=[pname '\'];
        fname=[fname ext];
    end
    
    load([pname fname]);
else
    fields=struct();
    fields.center=[0 0];
end

for i=1:length(fields)
    fields(i).fieldShift=-fields(i).center;
end

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
psfRange=10;%psf.range;
descr=psf.descr;

%target dx for making the psf. We will average later.
dx_t=.002;
dx_f=round(dx/dx_t);

%psf changes definition immediately below
npsf=round(psfRange./dx);
psf=zeros(round(psfRange./dx));
% [xpsf ypsf]=meshgrid([-npsf:1:npsf],[-npsf:1:npsf]);
% xpsf=xpsf.*dx;
% ypsf=ypsf.*dx;
% rpsf2=xpsf.^2+ypsf.^2;
% psf=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2)+eta/(pi*beta^2).*exp(-rpsf2./beta.^2));

%this is the upper right quadrant
%We need to start from -.5 to make sur the first points correspond to the
%psf at r=0;
[xpsf ypsf]=meshgrid([-.5:1:(npsf-1)*dx_f-1],[-.5:1:(npsf-1)*dx_f-1]);
xpsf=xpsf.*dx/dx_f;
ypsf=ypsf.*dx/dx_f;
rpsf2=xpsf.^2+ypsf.^2;
psf=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2)+eta/(pi*beta^2).*exp(-rpsf2./beta.^2));

%properly average the psf quadrant we have just created.
psf = reshape(psf, size(psf,1), dx_f,size(psf,2)/dx_f);
psf = squeeze(mean(psf, 2));
psf=psf';
psf=reshape(psf, size(psf,1), dx_f,size(psf,2)/dx_f);
psf = squeeze(mean(psf, 2));
psf=psf';

%Now we need to make the full psf by assembling all of the quadrants.
psf_ur=psf;
psf_ul=fliplr(psf(:,2:end));
psf_tot=[psf_ul, psf_ur];
psf_l=flipud(psf_tot(2:end,:));
psf_tot=[psf_l;psf_tot];
psf=psf_tot;

%Zero pad to at least 10um x 10 um;
%pad in the x direction
xpad=size(polysbin,1)-size(psf,1);
if xpad>0
    psf=padarray(psf,[xpad/2,0],0,'both');
elseif xpad<0
    polysbin=padarray(polysbin,[-xpad/2,0],0,'both');
end

%pad in the y direction
ypad=size(polysbin,2)-size(psf,2);
if ypad>0
    psf=padarray(psf,[0,ypad/2],0,'both');    
elseif ypad<0
    polysbin=padarray(polysbin,[0,-ypad/2],0,'both');
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
doseNew=dd(padPoints+1:end-padPoints-1*addY,padPoints+1:end-padPoints-1*addX);
shape=ss(padPoints+1:end-padPoints-1,padPoints+1:end-padPoints-1);
mp=size(doseNew,1);
np=size(doseNew,2);

fprintf('done.\n')

fprintf('Fracturing...\n')


%This is needed to not count the dose of places that get zero or NaN dose.
try
    doseNew(doseNew==0)=NaN;
    doseNew(doseNew<0)=NaN;
end

figure(557); clf;
hist(doseNew(:));
xlabel(dose);
drawnow;

doseStore=doseNew;

dvals=config.dvals;
nlayers=length(dvals);
dvalsl=dvals-(dvals(2)-dvals(1));
dvalsr=dvals;
layer=[];
figSize=ceil(sqrt(nlayers));
dvalsAct=[];


%loop variables used here
% j: fields (layers in the input dxf file)
% i: dose values (individual layers in the dxf files)
% n: subfield index
% m: subfield index
% b: boundaries
% k: boundaries with holes

%TODO:
% Loop over polygons
% fracture them individually
% 

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
    
    subField(ar).shape=p;
    subpoly = inpolygon(XPold, YPold, p(:,1), p(:,2));
    shotMap=doseNew.*subpoly;
    try
        shotMap(shotMap==0)=NaN;
        shotMap(shotMap<0)=NaN;
    end
    
    if fracture
        %Not really implemented yet.
        [b d]=fracturePoly(XP,YP,p,shotMap,dvals)
    else
        subField(ar).boundaries={p};
        subField(ar).dose=squeeze(nanmean(shotMap(:)));
        [mv,i]=min(abs(dvals-subField(ar).dose));
        subField(ar).doseColor=i;
        subField(ar).layer=layerNum(ar);
        subField(ar).field=1;
    end
    
    dvalsAct=dvals;

end

% for i=1:length(dvals)
%     fprintf('layer %d...\n',i);
%     
%     %Compute shot map for each layer
%     if i==1
%         layer(i).shotMap=doseNew<dvalsr(i).*shape;
%         layer(i).dose=dvals(i);
%         
%     elseif i==length(dvals)
%         layer(i).shotMap=doseNew>dvalsl(i).*shape;
%         layer(i).dose=dvals(i);
%         
%     else
%         layer(i).shotMap=(doseNew>dvalsl(i)).*(doseNew<dvalsr(i));
%         layer(i).dose=(dvalsl(i)+dvalsr(i))/2;
%     end
%     
%     %Compute the actual mean dose in the layer.
%     dd=doseNew;
%     dd(isnan(dd))=0;
%     doseSum= sum(sum(layer(i).shotMap));
%     if doseSum>0
%         layer(i).meanDose=sum(sum(layer(i).shotMap.*dd))/sum(sum(layer(i).shotMap));
%         dvalsAct(i)=layer(i).meanDose;
%     else
%         layer(i).meanDose=dvals(i);
%         dvalsAct(i)=dvals(i);
%     end
%     
%     %break into subfields and find boundaries. This is necessary because
%     %designCAD can only handle polygons with less than ~200 points.
%     
%     subfieldSize = config.subfieldSize;
%     
%     layer(i).boundaries={};
%     count=1;
%     
%     decrease = 0;
%     disp = 0;
%     m=1;
%     n=1;
%     xsubfields=ceil(mp/subfieldSize);
%     ysubfields=ceil(np/subfieldSize);
%     decrease_sub = 1;
%     while (n <= ysubfields)    %change to xsubfields for horizontal writing
%         %decrease_sub=1;
%         m = 1;
%         %decrease_sub = 1;
%         while (m <= xsubfields) %change to ysubfields for horizontal writing
%             %EJC: add decrease_sub factor for halving sub field size when polygons are large
%             if decrease
%                 subfieldSize = round(subfieldSize/decrease_sub);
%                 decrease = 0;
%                 disp = 0;
%             end
%             xsubfields=ceil(mp/subfieldSize);
%             ysubfields=ceil(np/subfieldSize);
%             if ~disp
%                 display(['trying subfield size of ' num2str(subfieldSize) '. There are a total of ' num2str(xsubfields*ysubfields) ' subfields...']);
%                 disp = 1;
%             end
%             proceed = 0;
%             
%             
%             subfield=zeros(mp,np);
%             %display(['Trying subfield size: ' num2str(mp/decrease_sub) 'x' num2str(np/decrease_sub)]);
%             if (m-1)*subfieldSize+1<min(m*subfieldSize,mp)
%                 xinds=(m-1)*subfieldSize+1:1:min(m*subfieldSize,mp);
%             else
%                 xinds=(m-1)*subfieldSize+1:1:min(m*subfieldSize,mp);
%             end
%             if (n-1)*subfieldSize+1 < min(n*subfieldSize,np)
%                 yinds=(n-1)*subfieldSize+1:1:min(n*subfieldSize,np);
%             else
%                 yinds=(n-1)*subfieldSize+1:1:min(n*subfieldSize,np);
%             end
%             
%             xstart=xinds(1);
%             ystart=yinds(1);
%             
%             %double the size of each shot map to avoid single pixel
%             %features. No longer used. It was an attemp to avoid
%             %single-pixel features.
%             %             xinds=reshape([xinds;xinds],[1 2*length(xinds)]);
%             %             yinds=reshape([yinds;yinds],[1 2*length(yinds)]);
%             
%             subdata=layer(i).shotMap(xinds,yinds);
%             
%             %Now "smear" out the shot map by one pixel in each direction. This makes sure that
%             %subfield boundaries touch each other.
%             sdll=padarray(subdata,[1,1],0,'pre');
%             sdur=padarray(subdata,[1,1],0,'post');
%             sdul=padarray(padarray(subdata,[1,0],'pre'),[0,1],'post');
%             sdlr=padarray(padarray(subdata,[1,0],'post'),[0,1],'pre');
%             
%             sd=sdll+sdur+sdul+sdlr;
%             sd(sd>0)=1;
%             subdata=sd;
%             
%             if (sum(subdata(:)))>0
%                 
%                 [B,L,nn,A]=bwboundaries(subdata);
%                 
%                 if debug
%                     figure(777); clf; hold on
%                     fprintf('Layer %d subfield (%d,%d) \n',i,m,n);
%                     for j=1:length(B)
%                         try
%                             bb=B{j};
%                             subplot(1,2,1); hold on;
%                             plot(bb(:,2),bb(:,1))
%                             axis([0 subfieldSize 0 subfieldSize]);
%                             subplot(1,2,2)
%                             imagesc(subdata); set(gca,'YDir','norm')
%                             axis([0 subfieldSize 0 subfieldSize]);
%                         end
%                         pause(1)
%                         
%                     end
%                     drawnow;
%                     pause(1);
%                 end
%                 
%                 if ~isempty(B)
%                     for b=1:length(B)
%                         
%                         %Find any holes and fix them by adding them
%                         %appropriately to enclosing boundaries
%                         enclosing_boundaries=find(A(b,:));
%                         for k=1:length(enclosing_boundaries)
%                             b1=B{enclosing_boundaries(k)};% the enclosing boundary
%                             b2=B{b}; %the hole
%                             
%                             %Only keep the hole if it has >0 area.
%                             %Also, only keep the hole if there are no
%                             %overlapping lines in the hole. There should be only one
%                             %pair of matching vertices per polygon.
%                             if polyarea(b2(:,1),b2(:,2))>0 && (size(b2,1)-size(unique(b2,'rows'),1)==1)
%                                 b1=[b1; b2; b1(end,:)]; %go from the enclosing boundary to the hole and back to the enclosing boundary
%                             end
%                             B{b}=[]; %get rid of the hole, since it is now part of the enclosing boundary
%                             B{enclosing_boundaries(k)}=b1; %add the boundary back to B.
%                         end
%                     end
%                     
%                     %Add boundaries to layer
%                     for b=1:length(B)
%                         if ~isempty(B{b}) && polyarea(B{b}(:,1),B{b}(:,2))>0
%                             
%                             %remove unnecessary vertices
%                             B{b}=simplify_polygon(B{b});
%                             
%                             %check for large polygons
%                             if size(B{b},1)>200%200
%                                 fprintf('Large boundaries in layer %d. Halving subfield size and retrying... \n',i);
%                                 
%                                 %If large boundaries, make the subfields
%                                 %smaller and restart;
%                                 decrease = 1;
%                                 decrease_sub = decrease_sub + 1;
%                                 proceed = 0;
%                                 disp = 0;
%                                 m = 1;
%                                 n = 1;
%                                 count = 1;
%                                 layer(i).boundaries = {};
%                                 %anybad = 1;
%                             elseif ~decrease
%                                 decrease_sub = 1;
%                                 proceed = 1;
%                             end
%                             
%                             if proceed
%                                 %divide by two because we doubled the size of the shot map
%                                 %subtract 1/2 because we smeared out the
%                                 %shot map by 1/2 of an original pixel
%                                 %Subtract 1/2 again because of the way we
%                                 %doubled the size of the matrix.
%                                 %layer(i).boundaries(count)={B{b}/2-1/2-1/2+repmat([xstart,ystart],[size(B{b},1),1])};
%                                 
%                                 %For use with undoubled matrices. Subtract
%                                 %1 because of the way we did the smearing
%                                 layer(i).boundaries(count)={B{b}-1+repmat([xstart,ystart],[size(B{b},1),1])};
%                                 
%                                 subField(m,n).boundaries(end+1)={B{b}-1+repmat([xstart,ystart],[size(B{b},1),1])};
%                                 subField(m,n).layer(end+1)=i;
%                                 
%                                 
%                                 count=count+1;
%                             end
%                         end
%                         
%                         %Break out of looping over boundaries if there are
%                         %large boundaries and we need to restart
%                         if proceed==0
%                             break
%                         end
%                     end
%                 end
%             else
%                 proceed = 1;
%             end
%             if proceed
%                 m = m+1;
%             end
%         end
%         if proceed
%             n = n+1;
%         end
%     end
%     
%     
%     fprintf('done.\n')
%     
%     
% end
fprintf('Fracturing complete. \n')

%Inspect layers if needed
if debug
    
    i=13;
    figure(777); clf; hold on
    
    maxSize=0
    for b=1:length(layer(i).boundaries)
        bb=layer(i).boundaries{b};
        plot(bb(:,2),bb(:,1))
        
        ms=size(layer(i).boundaries{b},1);
        if ms> maxSize
            maxSize=ms;
        end
    end
    maxSize
    
    figure(778); clf; imagesc(layer(i).shotMap)
    
end

%save the final files

%Coordinates for the pixels in the shot map
xpwrite=xpold;%[xpold xpold(end)+dx]-dx/2; %we possibly added one point to the array
ypwrite=ypold; %[ypold ypold(end)+dx]-dx/2; %we possible added one point to the array

fprintf('Fracturing write fields...')
%Go through the subFields, and find out which polygons belong to which fields
%TODO: need to loop through subfields here.
%JMN: not working yet
for i=1:length(layer)
    inds=ones(1,length(layer(i).boundaries));
    
    for j=1:length(fields)
        
        polygons=struct();
        polygons(1)=[];
        
        if isempty(inds)
            break
        end
        
        fields(j).layer(i).boundaries=[];
        
        for b=1:length(layer(i).boundaries)
            
            if ~inds(b)
                continue
            end
            
            xq=xpwrite(layer(i).boundaries{b}(:,2));
            yq=ypwrite(layer(i).boundaries{b}(:,1));
            if config.field
                xv=fields(j).box(:,1);
                yv=fields(j).box(:,2);
                
                [in,on] = inpolygon(xq,yq,xv,yv);
                
                if any(in) || any(on)
                    inds(b)=0; %its already in a field now
                    fields(j).layer(i).boundaries=[fields(j).layer(i).boundaries layer(i).boundaries(b)];
                    
                    
                end
            else
                fields(j).layer(i).boundaries=layer(i).boundaries;
            end
            
        end
    end
end

fprintf('done.\n')

for j=1:length(fields)
    
    outputFileName=[pathname filename(1:end-4) '_' descr '_' num2str(j) '.dxf'];
    dc2FileName=[pathname filename(1:end-4) '_' descr '_' num2str(j) '.dc2'];

    fprintf('Exporting to %s...\n',outputFileName);
    
    fields(j).cadFile=[filename(1:end-4) '_' descr '_' num2str(j) '.dc2'];
    
    FID = dxf_open(outputFileName);
    
    %ctab={[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1] [1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1] [1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [0 1 1]  };
    
    %Jet-like colors from working run file
    ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0] };
    
%     figure(558); clf; hold on;
%     title('Boundaries');
    %
    %     figure(559); clf; hold on;
    %     title('Shot Map');
    %     smap=doseNew.*0;
    
    polygons=struct();
    polygons(1)=[];
    
    %Write in order of layers
%     for i=length(dvals):-1:1
%         fprintf('Writing layer %d...\n',i)
%         %FID=dxf_set(FID,'Color',ctab{i}.*255,'Layer',i); % EJC: ctab{i} to ctab{i}.*255 (3/8/2019)
%         FID=dxf_set(FID,'Color',ctab{i}./255,'Layer',i); %JMN back to ctab{i}.*255*1. This is no longer needed because urpec saves dc2 files.
% 
%         %figure(558);
%         
%         for b=1:length(fields(j).layer(i).boundaries)
%             bb=fields(j).layer(i).boundaries{b};
%             X=xpwrite(bb(:,2)); X=X'+fields(j).fieldShift(1);
%             Y=ypwrite(bb(:,1)); Y=Y'+fields(j).fieldShift(2);
%             Z=X.*0;
%             dxf_polyline(FID,X,Y,Z);
%             %plot(X,Y);
%             
%             polygons(end+1).p=[X Y];
%             polygons(end).color=ctab{i}; %JMN 2019/10/12
%             polygons(end).layer=i;
%             polygons(end).lineType=1;
%             
%         end
%         %         figure(559);
%         %         smap=smap+layer(i).shotMap;
%         %         imagesc(smap);
%     end
    
    %Write in order of subfields. This makes more sense. Change the outer
    %loop to change the direction of writing.
    for ar=1:length(subField)
          
            for b=1:length(subField(ar).boundaries)
                i=subField(ar).doseColor(b);
                FID=dxf_set(FID,'Color',ctab{i}./255,'Layer',1); %JMN back to ctab{i}.*255*1. This is no longer needed because urpec saves dc2 files.
                bb=subField(ar).boundaries{b};
                X=bb(:,1); X=X+fields(j).fieldShift(1);
                Y=bb(:,2); Y=Y+fields(j).fieldShift(2);
                Z=X.*0;
                dxf_polyline(FID,X,Y,Z);
                %plot(X,Y);
                
                polygons(end+1).p=[X Y];
                polygons(end).color=ctab{i}; %JMN 2019/10/12
                polygons(end).layer=subField(ar).layer;
                polygons(end).lineType=1;
                
            end           
        
    end

    
    dxf_close(FID);
    
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
    

    
end

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

