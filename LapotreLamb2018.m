function [K_mean, K_max, K_min] = LapotreLamb2018(D50)
%% calculate hydraulic conductivity from empirical fits to natural data
% from Lapotre & Lamb (2018) Geology supplemental material
% inputs:
%   D50 = median grain size (m)
% outputs:
%   K_max = upper bound on hydraulic permeability (m2)
%   K_min = lower bound on hydraulic permeability (m2)

% constants
g = 9.81;                       % gravitational acceleration (m/s2)
nu = 10^(-6);                   % kinematic viscosity (m2/s)

% calculate lower and upper bounds on natural permeabilities
% take average to be representative value (Shepherd, 1989)
K_min = 11.9*(nu/g)*D50^1.5;
K_max = 6695*(nu/g)*D50^1.85;
K_mean = mean([K_min; K_max]);


end