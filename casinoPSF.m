function [  psf] = casinoPSF(  )
%function [  ] = casinoPSF(  )
%   Extracts a point spread function from a dataset generated using Casino.
%   http://www.gel.usherbrooke.ca/casino/What.html
%
%   You should run your simulation and save the data. Open the datafile in
%   excel and then save it as a .xlsx file. 
%
%   This function assumes that the top layer of your simulation is labeled
%   'PMMA' and calculates home much energy is deposited in the PMMA as
%   a function of position. You should simulate something like 10000
%   electrons. The excel file should be approximately ~75MB.
%   
%   TODO:
%   Add support for loading Casino files without intermediate excel step.

%load file
fprintf('Loading excel data file. This can take a while...');
[filename pathname]=uigetfile('*.xlsx');
[num txt raw]=xlsread([pathname filename]);
fprintf('done.\n');

% Find the energy deposited in the pmma vs lateral position.
fprintf('Counting electrons...');
dr=2; %spacing in nanometers for the simulation
rvals=(1:dr:10000);
evals=rvals.*0;
count=0; %counts the number of collisions in the pmma that were included in the loop.
pc=0; %counts the number of times an electron was in the pmma.
for i=1:length(raw)
    if strmatch(raw{i,8},'PMMA')
        pc=pc+1;
        try
            r=(raw{i,1}^2+raw{i,2}^2)^(0.5);
            ind=round(r/dr)+1;
            de=raw{i-1,7}-raw{i,7}; %energy deposted is difference in energy between successive collisions
            if length(de)==1 %the file format is weird. Sometimes de is not in the expected format.
                evals(ind)=evals(ind)+de;
                count=count+1;
            end
        catch
            %fprintf('fail \n');
            %ind;
        end
    end
    
    if strmatch(raw{i,2},'Trajectory')
        nelec=raw{i,3};
    end
end
% pc
% count
fprintf('done.\n');
fprintf('The average energy deposited per electron is %f keV. \n',sum(evals)/nelec);

%
%figure(333); clf; loglog(rvals,evals);

%Fit the data on a log log scale.
fitfn=@(p,x) p(1)/(1+p(2)).*((1/(pi*p(3)^2)).*exp(-x.^2/p(3).^2)+p(2)/(pi*p(4)^2)*exp(-x.^2/p(4).^2));

beta=[10000 .7 5 2000];

% figure(333); hold on;
% loglog(rvals,fitfn(beta,rvals));

% Fit the psf
%psf=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2)+eta/(pi*beta^2).*exp(-rpsf2./beta.^2));
%p(1) is amplitude
%p(2) is eta
%p(3) is alpha
%p(4) is beta
fitfn=@(p,x) log(p(1)/(1+p(2)).*((1/(pi*p(3)^2)).*exp(-(exp(x).^2)./p(3)^2)+p(2)/(pi*p(4)^2)*exp(-(exp(x).^2)./p(4)^2)));

%Divide by R because we calculated evals effectively by integrating rings.
ee=evals./rvals;

%Detemrine eta and keep it fixed during the fit
forwardRange=100; %We assume all forward scattering happens in less than 100 nm.
inds=(1:1:forwardRange/dr);
f=sum(evals(inds))-sum(evals(inds(end)+inds)); %approximation of the forward scattered energy
r=sum(evals(inds(end)+1:end)); %approximation of the back-scattered energy
eta=r/f;

%uncomment to keep an even density of points in logspace
% ii=linspace(0,log10(rvals(end)),200);
% ii=10.^(ii);
% ii=round(ii/dr);
% logKeep=zeros(1,length(rvals));
% logKeep(ii)=1;

logKeep=ones(1,length(rvals));

y=log(ee);
reject=(y==-Inf);
keep=~reject;
keep=keep & logKeep;
y=y(keep);
x=log(rvals);
x=x(keep);
beta0=[20000 eta 2 1448];
beta=fitwrap('plinit plfit robust',x,y,beta0,fitfn,[1 0 1 1]);
xlabel('log(range [nm])');
ylabel('log(energy deposited per distance [kV/nm])');
title('Data and fit');

beta=real(beta);

psf=struct();

psf.eta=beta(2); %ratio of backscattered energy to forward scattered energy
psf.alpha=beta(3).*1e-3;  %forward scattering range, in units of microns
psf.beta=beta(4).*1e-3;  %backward scattering range, in units of microns
psf.range=5;
psf.x=x;
psf.y=y;
psf.beta0=beta0;
psf.fitfn=fitfn;
psf.mask=[1 0 1 1];
descr=input('Enter a description, e.g., Si30kV \n','s');
psf.descr=descr;

outputFileName=['PSF' descr '.mat'];
fprintf('Saving to %s',outputFileName);
save(outputFileName,'psf');

% Previously calculated PSF parameters

% %GaAs 30 kV with 200 nm pmma:
% eta=1.7
% alpha=14 nm
% beta=1.56 um
% 
% 
% % GaAs 30 kV with 400 nm pmma
% eta=1.5
% alpha=28 nm
% beta=1.74 um
% 
% % Si 30 kV with 200 nm pmma
% eta=1
% alpha=13 nm
% beta=3.97 um
% 
% % Si 30 kV with 400 nm pmma
% eta=0.9
% alpha=40 nm
% beta=4.19 um

end

