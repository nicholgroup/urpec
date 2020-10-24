pecdir=fileparts(which('urpec'));
curdir=pwd;
cd(pecdir);

%% Easy test, default, but the dxf file has duplicate points for each shape.
% There are no concave polygons in this pattern.
psf='PSFLiNbO3_30kV_200.mat';
filenameP='saw.dxf';
pathnameP='Examples\'

urpec_v3(struct('file',[pathnameP filenameP],...
    'psfFile',psf));


%% Fracturing, default
psf='PSFSi90SiO230kV200.mat';
filenameP='2D.dxf';
pathnameP='Examples\'

urpec_v3(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%% Regular, wrong layer names
psf='PSFGaAs30kV200.mat';
filenameP='P7.dxf';
pathnameP='Examples\'

urpec_v3(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%% Off center, wrong layer names
psf='PSFGaAs30kV200.mat';
filenameP='MM.dxf';
pathnameP='Examples\'

urpec_v3(struct('file',[pathnameP filenameP],...
    'psfFile',psf,...
    'outputDir','C:\Users\Nichol\Desktop\tmp'));


%% Small features, off center
psf='PSFSi30kV200.mat';
filenameP='dot.dxf';
pathnameP='Examples\'

urpec_v3(struct('file',[pathnameP filenameP],...
    'autoRes',false,...
    'dx',.001,...
    'dvals',linspace(1,2.0,15),...
    'padLen',.5,...
    'maxIter',4,...
    'fracNum',4,...
    'fracSize',4,...
    'psfFile',psf));
