function [e] = psf3(p,x)
%psf3 computes the density of scattered energy in electron beam lithography
%   
% This is a fit function. All distance argument are logarithms of actual values.
% See Proximity Effect in E-beam Lithography at
% 1 for more details

alpha=10^p(1); %alpha, forward scattering range
beta=10^p(2); %beta, reverse scattering range
gamma=10^p(3); %gamma, exponential tail range

et=p(4); %ratio of total back scattered energy to forward-scattered energy.
eta=et-p(5); %eta
nu=p(5); %gamma

r=10.^x; %distance from the center

e= log10(1/(pi*(1+eta+nu)).*((1/alpha^2).*exp(-(r.^2)./(alpha^2))+ (eta/beta^2).*exp(-(r.^2)./(beta^2))+(nu/(24*gamma^2)).*exp(-(r./gamma).^(1/2))));


end

