function [fD, Cf] = CWDarcyWeisbach(H, U, ks)
%% calculate the Darcy-Weisbach friction factor and coefficient of friction from the Colebrook-White approximation for transitional flow in an open channel
% inputs:
%   H = flow depth (m)
%   U = average flow velocity (m/s)
%   ks = roughness height (m)

nu = 1e-6;                      % kinematic viscosity of water (m2/s)
g = 9.81;                       % gravitational acceleration (m2/s)
Re = U * H / nu;                % fluid Reynolds number (unitless)

CW = @(x) 1/sqrt(x) + 2*log10(2.51/Re/sqrt(x) + ks/(3.7*H));
fD = fsolve(CW,0.1);
Cf = fD/8;

end