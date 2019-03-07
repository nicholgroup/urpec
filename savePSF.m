%% 
% This is a script to manually save data from casino.
% It is better to run casinoPSF;

%% Si 30 kV from Casino
fname='Si30kVPSF.mat'
psf=struct();
psf.eta=0.67; %ratio of backscattered energy to forward scattered energy
psf.alpha=0.01;  %forward scattering range, in units of microns
psf.beta=2.85;  %backward scattering range, in units of microns
psf.range=5;
psf.descr='Si30kV';

save(fname,'psf');

% eta=0.67; 
% alpha=.01;
% beta=2.85;
% psfRange=5;

%% GaAs 30 kV from Casino
fname='GaAs30kVPSF.mat'
psf=struct();
psf.eta=1.7; %ratio of backscattered energy to forward scattered energy
psf.alpha=0.014;  %forward scattering range, in units of microns
psf.beta=1.56;  %backward scattering range, in units of microns
psf.range=5;
psf.descr='GaAs30kV';

save(fname,'psf');

% %GaAs 30 kV, based on Casino
% eta=1.7; %ratio of backscattered energy to forward scattered energy
% alpha=.014; %forward scattering range, in units of microns
% beta=1.56; %backward scattering range, in units of microns
% psfRange=5;