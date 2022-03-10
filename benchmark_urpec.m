%Contains different excercises for urpec. After changes to urpec, it needs
%to pass these tests.

pecdir=fileparts(which('urpec'));
curdir=pwd;
cd(pecdir);


%% Terrible, non-convex polygon
filenameP='2d_v2.dxf';
pathnameP='Examples\';

fieldsFile=urpec_v4(struct('file',[pathnameP filenameP],...
    'autoRes',true,...
    'dx',1,...
    'psfFile','PSFSi90SiO250kV200.mat'));

%Display the pattern with doses
plotFields(fieldsFile);

%% Easy test, default, but the dxf file has duplicate points for each shape.
% There are no concave polygons in this pattern.
psf='PSFLiNbO3_30kV_200.mat';
filenameP='saw.dxf';
pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%% Lots of polygons, speed test.

psf='PSFLiNbO3_30kV_200.mat';
filenameP='f500M_IDTp1_REFp200_GJ_wo2_C.mat';
pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf));


%% Fracturing, default
psf='PSFSi90SiO250kV200.mat';
filenameP='2D_v2.dxf';
pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'dx',1,...
    'psfFile',psf));

%% Regular, wrong layer names
psf='PSFGaAs30kV200.mat';
filenameP='P7.dxf';
pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%% Off center, wrong layer names, very fast
psf='PSFGaAs30kV200.mat';
filenameP='MM.dxf';
pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf,...
    'autoRes',false,...
    'dvals',linspace(1,2.0,32),...
    'npgs',true,...
    'dx',.05));

%% Spiral tests
psf='PSFSiNb125_950_50.mat';
%filenameP='spiral_nice.dxf';
%filenameP='spiral_nice.mat';
filenameP='spiral_bad.mat';

pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf,...
    'autoRes',false,...
    'dx',.05));


%% Small features, off center. 
%This can easily generate lots of points, so be careful
psf='PSFSi30kV200.mat';
filenameP='dot.dxf';
pathnameP='Examples\';

urpec_v4(struct('file',[pathnameP filenameP],...
    'autoRes',false,...
    'dx',.005,...
    'dvals',linspace(1,2.0,32),...
    'padLen',.5,...
    'maxIter',4,...
    'fracNum',4,...
    'fracSize',4,...
    'psfFile',psf));
