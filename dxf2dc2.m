function [  ] = dxf2dc2( config )
% function [  ] = dxf2dc2( config )
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
%   file: datafile for processing. This can either be a .dxf file or a .mat
%
% By:
% John Nichol jnich10@ur.rochester.edu
%


if ~exist('config','var')
    config=struct();
end

config=def(config,'file',[]); 

ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0] };

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
    
if strmatch(ext,'.dxf')
    [lwpolylines,lwpolylayers]=dxf2coord_20(pathname,filename);
    %For later use in breaking into fields
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
    layerNum=[d.polygons.layer];
    lwpolylines=[ones(size(d.polygons(1).p,1),1).*1 d.polygons(1).p];
    for i=2:length(d.polygons)
        lwpolylines=[lwpolylines; [ones(size(d.polygons(i).p,1),1).*i d.polygons(i).p]];
    end
end

fprintf('Converting...');

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

%Make a final polygon array for saving
polygons=struct();
polygons(1)=[];

%Write in order of objects.
for ar=1:length(objects)
    
    X=objects{ar}(:,2);
    Y=objects{ar}(:,3);
    
    polygons(end+1).p=[X Y];
    polygons(end).color=ctab{layerNum(ar)}; 
    polygons(end).layer=layerNum(ar);
    polygons(end).lineType=1;
    
end

dc2FileName=[pathname filename(1:end-4) '.dc2'];

dc2write(polygons,dc2FileName);

fprintf('done. \n');


end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
end

