%% Setup
files={'C:\Users\Nichol\Box Sync\Nichol Group\Nichol\fab\SiGe_Qubit3_L3_sm.dxf',...
    'C:\Users\Nichol\Box Sync\Nichol Group\Nichol\fab\SiGe_Qubit3_L3_med.dxf',...
    'C:\Users\Nichol\Box Sync\Nichol Group\Nichol\fab\SiGe_Qubit3_L3_lg.dxf'};

psfFile='C:\Users\Nichol\Documents\GitHub\urpec\PSFSi30kV200.mat';

config=struct();
config.file=files{1};
config.psfFile=psfFile;
config.dvals=[1:.1:1.9];

%% run

for i=1:length(files)
    config.file=files{i};
    urpec(config);
end
    
