function [Usub1, Usub2, P] = Lamb2017WRR_limited(H, S, D50, K, Lambda)
%% calculate the average subsurface flow velocity for flow through a porous medium with an applied surface shear stress
% inputs:
%   H = flow depth (m)
%   S = open channel slope (m/m)
%   D50 = median grain size (m)
%   K = hydraulic permeability (m2)
%   Lambda = subsurface porosity (-)
% outputs:
%   Usub1 = mean near-surface flow velocity (m/s)
%   Usub2 = mean deep flow velocity (m/s)
%   P = exchange layer depth (m)

% hydraulic constants
rho = 1000;                         % density of water (kg/m3)
mu = 0.0015;                        % dynamic viscosity (Pa*s)
g = 9.81;                           % gravitational acceleration (m/s2)

k = rho*g*K/mu;                     % hydraulic conductivity (-)   
F = 5*10^(-3);                      % Forcheimer coefficient (-)
F_star = F * Lambda / g / sqrt(K);         % (s2/m2)  

% calculate exchange layer depth P (m)
CD = 1;                             % Nepf et al. (2012), bed drag coefficient
a = 3*(1 - Lambda)/(2*D50);         % Ghisalberti (2009), xs area (m2)
P = 1/(CD * a * 3);                 % Ghisalberti (2009)

% calculate mean flow velocity for upper layer of sediment (m/s)
C1 = 2;             % linear velocity profile with depth in boundary layer
Usub1 = -1/(2*C1*F_star*k) + 0.5*sqrt((k*C1*F_star)^(-2) + 4*S*(1 + H/P)/(C1*F_star));

% calculate mean flow velocity for lower layer of sediment (m/s)
% Usub2 = -1/(2*C1*F_star*k) + 0.5*sqrt((k*C1*F_star)^(-2) + 4*Grad/(C1*F_star));
Usub2 = 0;

end