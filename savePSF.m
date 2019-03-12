%% 
% This is a script to manually save data from casino.
% It is better to run casinoPSF;

%% Si 30 kV 200nm PMMA from Casino
fname='PSFSi30kV200.mat'
psf=struct();
psf.eta=0.67; %ratio of backscattered energy to forward scattered energy
psf.alpha=0.01;  %forward scattering range, in units of microns
psf.beta=2.85;  %backward scattering range, in units of microns
psf.range=5;
psf.resist=200; %resist thickness, nm
psf.descr='Si30kV200';

save(fname,'psf');

% eta=0.67; 
% alpha=.01;
% beta=2.85;
% psfRange=5;

%% GaAs 30 kV 200 nm PMMA from Casino
fname='PSFGaAs30kV200.mat'
psf=struct();
psf.eta=1.7; %ratio of backscattered energy to forward scattered energy
psf.alpha=0.014;  %forward scattering range, in units of microns
psf.beta=1.56;  %backward scattering range, in units of microns
psf.range=5;
psf.resist=200; %resist thickness, nm
psf.descr='GaAs30kV200';

save(fname,'psf');

% %GaAs 30 kV, based on Casino
% eta=1.7; %ratio of backscattered energy to forward scattered energy
% alpha=.014; %forward scattering range, in units of microns
% beta=1.56; %backward scattering range, in units of microns
% psfRange=5;

%% GaAs 30 kV 400 nm PMMA from Casino
fname='PSFGaAs30kV400.mat'
psf=struct();
psf.eta=1.5; %ratio of backscattered energy to forward scattered energy
psf.alpha=0.028;  %forward scattering range, in units of microns
psf.beta=1.74;  %backward scattering range, in units of microns
psf.range=5;
psf.resist=400; %resist thickness, nm
psf.descr='GaAs30kV400';

save(fname,'psf');

% % GaAs 30 kV with 400 nm pmma
% eta=1.5
% alpha=28 nm
% beta=1.74 um
% 