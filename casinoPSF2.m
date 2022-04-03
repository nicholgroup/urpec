function [psf] = casinoPSF2(config)
%casinoPSF2 extracts a point spread function from a Casino dataset. 
%
%   Extracts a point spread function from a dataset generated using Casino.
%   http://www.gel.usherbrooke.ca/casino/What.html
%
%   Run the casino simulation and save the data as a .dat file
%
%   This function assumes that the top layer of your simulation has the
%   string 'MMA' in it. labeled and calculates the distribution of
%   electron collisions as af unction of radial distance. 
%   a function of position. You should simulate something like 2,000-10,000
%   electrons. The dat file will be about 150-600 MB.

if ~exist('config','var')
    config=struct();
end

config=def(config,'beta0',[1 4  1.1 0 .1]); %parameters are [log10(alpha), log10(beta), log10(gamma), total backscattered ratio, nu]
config=def(config,'file',[]); 

[filename pathname]=uigetfile({'*.dat';'*.mat'});
[~,~,ext] = fileparts(fullfile(pathname,filename));

%analyze a new simulation
if strmatch(ext,'.dat')
    fID=fopen([pathname filename]);
    count_elec=1;
end

%in case we want to reanalyze data.
if strmatch(ext,'.mat')
    count_elec=0;
    d=load([pathname filename]);
    rvals=d.psf.rvals;
    evals=d.psf.evals;
end

dr=1; %spacing in nanometers for the simulation
if count_elec
    % Find the energy deposited in the pmma vs lateral position.
    fprintf('Counting electrons\n');
    rvals=(1:dr:50000);
    evals=rvals.*0;
    nvals=rvals.*0;
    count=0; %counts the number of collisions in the pmma that were included in the loop.
    xv=[];
    yv=[];
    maxE=0;
    lastLine=[];
    
    cont=1;
    tic
    while cont
        A=fgetl(fID);
        
        if A==-1
            cont=0;
            break
        end
        
        skip=isempty(A);
        if skip
            continue
        end
        
        [aa,n]=sscanf(A,'%f %f %f %f %f %f %f %s');
        
        if n==8
            try
                de=lastLine(7)-aa(7);
            catch
                de=0;
            end
            
            if ~isempty(regexp(A,'MMA')) %Only look at collisions in PMMA or MMA
                xv=[xv aa(1)];
                yv=[yv aa(2)];
                r=(aa(1)^2+aa(2)^2)^(0.5);
                ind=round(r/dr)+1;
                
                if ind<length(rvals)
                    nvals(ind)=nvals(ind)+1; %just count collisions
                    evals(ind)=evals(ind)+de; %count energy
                end
                
                count=count+1;
            end
            
        end
        
        lastLine=aa;
        
        if count>0 & ~mod(count,1000)
            fprintf('%d collisions \n',count);
        end
        
    end
    toc
    fclose(fID);
end

r2=pi.*rvals.^2; %circle area
aa=r2-[0 r2(1:end-1)];  %right area

evals=evals./sum(evals); %normalize
edensity=evals./aa; %convert to density, with units of 1/nm^2

figure(444); clf; loglog(rvals,edensity,'o');
fprintf('Click on forward scattering boundary \n');
f=ginput(1);

%Determine the ratio of backward to forward scattering and keep it fixed during the fit
forwardRange=f(1); 
inds=(1:1:forwardRange/dr);
f=sum(evals(inds));%;-sum(evals(inds(end)+inds)); %approximation of the forward scattered energy, optionally trying to subtract a back-scattered background
r=sum(evals(inds(end)+1:end)); %approximation of the back-scattered energy
et=r/f; %ratio of backscattered energy to forward energy

fitfn=@(p,x) psf3(p,x);
beta0=config.beta0;
beta0(1)=log10(forwardRange/3); %parameters are [log10(alpha) log10(beta) log10(gamma) total backscattered ratio nu]
beta0(4)=et;  
mask=[1 1 1 0 1];

logKeep=ones(1,length(rvals));
logKeep(1:2)=0;

%Get rid of any points with no counts
y=log10(edensity);
reject=(y==-Inf);
keep=~reject;
keep=keep & logKeep;
y=y(keep);
x=log10(rvals);
x=x(keep);

beta=fitwrap('plinit plfit ',x,y,beta0,fitfn,mask);
beta=real(beta);

a=10^beta(1);
b=10^beta(2);
g=10^beta(3);
nu=beta(5);
eta=et-nu;

xlabel('log10(range [nm])');
ylabel('log10(energy deposited [keV])');
title(sprintf('${\\alpha}$=%3.3f nm, ${\\beta}$=%3.3f nm,  ${\\eta}$=%3.3f, ${\\gamma}$=%3.3f nm, ${\\nu}$=%3.3f',a,b,eta,g,nu),'interpreter','Latex');
box on;

psf=struct();

psf.params=beta;
psf.fitfn=fitfn;
psf.beta0=beta0;
psf.mask=mask;
psf.x=x;
psf.y=y;
psf.rvals=rvals;
psf.evals=evals;
psf.edensity=edensity;

descr=input('Enter a description, e.g., Si30kV_500 \n','s');
psf.descr=descr;
psf.version=2; 

outputFileName=['PSF' descr '.mat'];
fprintf('Saving to %s',outputFileName);
save(outputFileName,'psf');

end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
end
