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

fieldsFile=urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%display the fracture pattern
plotFields(fieldsFile);


%% Step 3: Make a default run file for NPGS
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
fieldsFile=urpec_v4(struct('file',[pathnameP filenameP],...
    'autoRes',false,...
    'dx',.001,...
    'dvals',linspace(1,2.5,32),...
    'maxIter',4,...
    'fracNum',3,...
    'fracSize',4,...
    'psfFile',psf));

% fieldsFile=urpec_v3(struct('file',[pathnameP filenameP],...
%     'dvals',linspace(1,2.5,15),...
%     'psfFile',psf));

plotFields(fieldsFile);