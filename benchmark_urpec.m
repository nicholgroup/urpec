%Contains different excercises for urpec. After changes to urpec, it needs
%to pass these tests.

pecdir=fileparts(which('urpec'));
curdir=pwd;
cd(pecdir);
T=[];

%% urpec logo
psf='PSFGaAs30kV200.mat';
filenameP='logo.mat';
pathnameP='Examples\';

[fieldsFile,T(1)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'convexify',false,...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);


%% Terrible, non-convex polygon
filenameP='2d_v2.dxf';
pathnameP='Examples\';

[fieldsFile,T(2)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'autoRes',true,...
    'convexify',false,...
    'dx',.15,...
    'psfFile','PSFSi90SiO250kV200.mat'));

%Display the pattern with doses
plotFields(fieldsFile);

%% Terrible, non-convex polygon with convexify option
filenameP='2d_v2.dxf';
pathnameP='Examples\';

[fieldsFile,T(3)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'autoRes',true,...
    'dx',.15,...
    'convexify',true,...
    'psfFile','PSFSi90SiO250kV200.mat'));

%Display the pattern with doses
plotFields(fieldsFile);

%% Easy test, default, but the dxf file has duplicate points for each shape.
% There are no concave polygons in this pattern.
psf='PSFLiNbO3_30kV_200.mat';
filenameP='saw.dxf';
pathnameP='Examples\';

[fieldsFile,T(4)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'convexify',true,...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);

%% Lots of polygons, speed test.

psf='PSFLiNbO3_30kV_200.mat';
filenameP='f500M_IDTp1_REFp200_GJ_wo2_C.mat';
pathnameP='Examples\';

[fieldsFile,T(5)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);

%% Fracturing, default
psf='PSFSi90SiO250kV200.mat';
filenameP='2D.dxf';
pathnameP='Examples\';

[fieldsFile,T(6)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'dx',.15,...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);

%% Regular, wrong layer names
psf='PSFGaAs30kV200.mat';
filenameP='P7.dxf';
pathnameP='Examples\';

[fieldsFile,T(7)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);

%% Off center, wrong layer names, very fast
psf='PSFGaAs30kV200.mat';
filenameP='MM.dxf';
pathnameP='Examples\';

[fieldsFile,T(8)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf,...
    'autoRes',false,...
    'dvals',linspace(1,2.0,32),...
    'npgs',true,...
    'dx',.05));
%Display the pattern with doses
plotFields(fieldsFile);

%% Spiral tests
psf='PSFSiNb125_950_50.mat';
%filenameP='spiral_nice_simple.dxf';
filenameP='spiral_nice.mat';
%filenameP='spiral_bad.mat';

pathnameP='Examples\';

[fieldsFile,T(9)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'psfFile',psf,...
    'autoRes',true,...
    'dx',.05));
%Display the pattern with doses
%plotFields(fieldsFile);

%% Small features, off center. 
%This can easily generate lots of points, so be careful
psf='PSFSi30kV200.mat';
filenameP='dot.dxf';
pathnameP='Examples\';

[fieldsFile,T(10)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'autoRes',false,...
    'dx',.005,...
    'dvals',linspace(1,2.0,32),...
    'padLen',.5,...
    'maxIter',4,...
    'fracNum',10,...
    'fracSize',2,...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);
%% Dxf has both 2D and 3D polylines present
psf='PSFGaAs30kV200.mat';
filenameP='3D_2D.dxf';
pathnameP='Examples\';

[fieldsFile,T(11)]=urpec_v4(struct('file',[pathnameP filenameP],...
    'convexify',false,...
    'psfFile',psf));

%Display the pattern with doses
plotFields(fieldsFile);
%%
round(T')
fprintf('Passed all tests \n');

