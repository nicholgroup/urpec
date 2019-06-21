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

%% Make a run file

config=struct();

config.dtc='300';
config.sm_aperture=30;
config.lg_aperture=60;
config.write_mag_sm='1100';

urpec_writeJob(config);
        

%% Make a run file for the SAW. 
    
