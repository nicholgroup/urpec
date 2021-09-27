%% README
% Run steps 1-3 to run urpec and make a run file. 

% The layer scheme for your cad file is as follows. 
% The names for all layers should be numbers.
% Layers 1 and 2 of the input file will
% both be output to layer 1 of the output file. Layer 1 will not be
% fractured, and layer 2 will be fractured. Layers 3 and 4 of the input
% file will be output to layer 2 of the output filed, etc. If the polygons
% are not fractured, they are written with an average dose. 
%
% Any layer names that do not respect this format will be fractured and
% output to layer 1 of the dc2 file.
%% Step 1: Choose pattern file and psf file
% This works best if your pattern file is already in the
% NPGS project directory.

[filenameP pathnameP]=uigetfile({'*.dxf';'*.mat'},'Select pattern file.');

[psf ]=uigetfile('PSF*.*','Select PSF file.');

%% Step 2: Run urpec with default settings

urpec_v3(struct('file',[pathnameP filenameP],...
    'psfFile',psf));


%% Step 3: Make a default run file
% Before you run this cell, make sure the field files ('...fields.mat'), dose files ('.txt'), 
% and the pattern files ('.dc2') are in the proper NPGS project directory. 
% If your pattern file was already in this directory before pec, everything
% should be ok.
%
% It should be possible to run the run file right away, but sometimes you
% may have to resave the .dc2 file for NPGS, and the reload and resave the
% run file.

%Magnification
mag=1500;

%aperture
ap=10; %can be 7,10,30,120

%input the current for the aperture.
current=40;

%dose
dtc=[400];

%L-L and C-C spacing, in Angstroms.
spacing=[50.743];

%Initial move
moves={[300,300]};

%pattern files
[filename pathname]=uigetfile('*fields.mat');
files1={filename};

%write layers
layers=[1];

dir=pwd;
cd(pathname);

% Create the structure of the run file here. 
entities=[];
entities.val={};

entities(1).type='header';
entities(1).val={};
entities(1).dir=pathname;

entities(2).type='magComment';
entities(2).val={mag};

pos=[0 0];

%You can provide array or structs for most of the parameters above, and you
%can define a run file that is many entities long.
for i=1:length(files1)
    
    d=load(files1{i});
    fields=d.fields;
    
    startInd=2+(i-1)*3+1;
    
    entities(end+1).type='move';
    entities(end).val={moves{i}(1),moves{i}(2)};
    
    %specity important write parameters here.
    entities(end+1).type='write';
    entities(end).val={fields(1).cadFile(1:end-4)};
    entities(end).dtc = num2str(dtc(i));
    entities(end).mag = num2str(mag);
    entities(end).aperture = ap;
    entities(end).current=current;
    entities(end).cadFile=fields(1).cadFile;
    entities(end).doseFile=fields(1).doseFile;
    entities(end).spacing={num2str(spacing(i)) ,num2str(spacing(i))};
    entities(end).layer=layers(i);
    
    pos=pos+moves{i};
   
end

entities(end+1).type='move';
entities(end).val={-pos(1),-pos(2)};

for i=1:length(entities)
    entities(i).dir=pathname;
end

entities(1).val={length(entities)-1};

%dir = pwd;
urpec_makeRunFile_v3(entities);

cd(dir);


%% Run urpec with advanced settings
% Here you can specify different parameters.
urpec_v3(struct('file',[pathnameP filenameP],...
    'autoRes',false,...
    'dx',.001,...
    'dvals',linspace(1,2.0,15),...
    'maxIter',4,...
    'fracNum',3,...
    'fracSize',4,...
    'psfFile',psf));


%% ########## Deprecated below ########## 

%% Set up for basic run
files={'C:\Users\Nichol\Box Sync\Nichol Group\Nichol\fab\SiGe_Qubit3_L3_sm.dxf',...
    'C:\Users\Nichol\Box Sync\Nichol Group\Nichol\fab\SiGe_Qubit3_L3_med.dxf',...
    'C:\Users\Nichol\Box Sync\Nichol Group\Nichol\fab\SiGe_Qubit3_L3_lg.dxf'};

psfFile='C:\Users\Nichol\Documents\GitHub\urpec\PSFSi30kV200.mat';

config=struct();
config.file=files{1};
config.psfFile=psfFile;
config.dvals=[1:.1:1.9];

%% Run for mulitple files

for i=1:length(files)
    config.file=files{i};
    urpec(config);
end

%% Run urpec for a single file, do not break into layers.
urpec(struct('file','C:\Users\Nichol\Box Sync\Nichol Group\Nichol\matlab\saw\2019_06_23_1G_p\Test.dxf','psfFile','PSFGaAs30kV200.mat'));

%% Run urpec for a file with separate layers
%Note: you need to define the layer shift appropriately for your pattern
curDir=pwd;
cd('C:\Users\Nichol\Documents\GitHub\urpec');

layerShift=[];
for i=1:length(layers)
    layerShift=[layerShift; -layers(i).center];
end

urpec(struct('file','C:\Users\Nichol\Box Sync\Nichol Group\Nichol\matlab\saw\saw.dxf',...
    'layer',true,...
    'layerShift',layerShift,...
    'psfFile','PSFGaAs30kV200.mat'));

cd(curDir);

%% Make a run file based on entities

%manually make a field struct, if urpec didn't already make it

% fields=struct();
% fields(1).cadFile='sawp_gaas30kv200_1.dc2';
% fields(1).doseFile='sawp_gaas30kv200_1.txt';
% fields(1).layerShift=[-133.534,0];
% 
% fields(2).cadFile='sawp_gaas30kv200_2.dc2';
% fields(2).doseFile='sawp_gaas30kv200_2.txt';
% fields(2).layerShift=[-276.734,0];
% 
% fields(3).cadFile='sawp_gaas30kv200_3.dc2';
% fields(3).doseFile='sawp_gaas30kv200_3.txt';
% fields(3).layerShift=[133.534,0];
% 
% fields(4).cadFile='sawp_gaas30kv200_4.dc2';
% fields(4).doseFile='sawp_gaas30kv200_4.txt';
% fields(4).layerShift=[276.734,0];
% 
% fields(5).cadFile='sawp_gaas30kv200_5.dc2';
% fields(5).doseFile='sawp_gaas30kv200_5.txt';
% fields(5).layerShift=[0,.122];
% 
% %save('sawp_gaas30kv200_fields.mat','fields');

% load the fields struct that urpect outputs. 
[filename pathname]=uigetfile('*fields.mat');
load(fullfile(pathname,filename));

p=[];
for i=1:length(fields)
    fields(i).position=-fields(i).layerShift;
    p=[p;-fields(i).layerShift];
end

[~,layerOrder]=sort(p(:,1));
positions=p(layerOrder,:);
moves=diff(positions);
%assumes you started from the pattern center
moves=[positions(1,:);moves];

initMove=[350,350];

% Create the structure of the run file here. 

entities=[];
entities.val={};

entities(1).type='header';
entities(1).val={};
entities(1).dir=pathname;

entities(2).type='magComment';
entities(2).val={650};

entities(3).type='move';
entities(3).val={initMove(1),initMove(2)};

%setions to repeat start here

entities(4).type='move';
entities(4).val={-400,400};

entities(5).type='waitComment';
entities(5).val={};

entities(6).type='align';
entities(6).val={'align','A1'};
entities(6).mag = '650';
entities(6).aperture = 30;
entities(6).current=450; 
entities(6).cadFile='align.dc2'; 

entities(7).type='write';
entities(7).val={'saw_gaas30kv200_1'};
entities(7).dtc = '230';
entities(7).mag = '650';
entities(7).aperture = 30;
entities(7).current=450; 
entities(7).cadFile='saw_gaas30kv200_1.dc2'; 
entities(7).doseFile='saw_gaas30kv200_1.txt'; 
entities(7).spacing={'97.6093' ,'97.6093'}; 

%sections to repeat end here

startInd=4; endInd=7; len=endInd-startInd+1;

%loop over the write fields, and populate write and align blocks
for i=1:length(fields)
    entities(4+(i-1)*len)=entities(4);
    entities(4+(i-1)*len).val={moves(i,1),moves(i,2)};
    
    entities(5+(i-1)*len)=entities(5);
    
    entities(6+(i-1)*len)=entities(6);
    entities(6+(i-1)*len).dir=pathname;

    
    entities(7+(i-1)*len)=entities(7);
    entities(7+(i-1)*len).val={fields(layerOrder(i)).cadFile(1:end-4)};
    entities(7+(i-1)*len).cadFile=fields(layerOrder(i)).cadFile;
    entities(7+(i-1)*len).doseFile=fields(layerOrder(i)).doseFile;
   
    entities(7+(i-1)*len).dir=pathname;
 
end

entities(end+1).type='move';
lastMove=-(sum(moves)+initMove);
entities(end).val={lastMove(1),LastMove(2)};

entities(1).val={length(entities)-1};

urpec_makeRunFile(entities);

%% Make a run file with urpec_writeJob

config=struct;
config.dtc = '230';
config.sm_aperture = 7.5;
config.write_mag_sm = '1500';

urpec_writeJob(config);
