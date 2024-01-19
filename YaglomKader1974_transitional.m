function [qw, Ch] = YaglomKader1974_transitional(U, H, Hr, Cf, Twater, Tbank)
%% calculate heat flux from fully turbulent fluid into rough boundary
% inputs:
%   U = average open channel flow velocity (m/s)
%   H = flow depth (m)
%   Hr = height of roughness elements into flow (m)
%   Cf = coefficient of friction (-)
%   Twater = temperature of channel water (K)
%   Tbank = temperature of bank surface (K)
% output:
%   qw = heat flux from flow to bank (J/m2/s)

% empirical constants from YK74
alpha = 2.12;
beta1 = 0.05;
C = 9.5;
b1 = 0.55;
b2 = 1/11;

% general empirical constants
g = 9.81;               % gravitational acceleration (m/s2)
Pr = 10;                % Prandtl number (-)
nu = 10^(-6);           % kinematic viscosity (m2/s)
cp = 4184;              % specific heat of water (J/kg/K)
rho = 1000;             % density of water (kg/m3)

n1 = Hr / H;                        % relative roughness height (-)
u_star = U*sqrt(Cf);                % open channel shear velocity (m/s)
h = Hr * u_star / nu;               % non-dimensional roughness height (-)
beta_s = 12.5*Pr^(2/3) - 6;
beta_r = sqrt(h)*(b1*Pr^(2/3) - b2) + C;

% calculate transitional flow Stanton number
if h > 100
    beta_t = beta_r;
else
    beta_t = (h/100)*beta_r + (1 - h/100)*beta_s;
end
Ch = sqrt(Cf) / (-alpha*log(n1) + beta_t + beta1);
% calculate heat flux to wall using equation 14
qw = cp * rho * U * Ch * (Twater - Tbank);


end