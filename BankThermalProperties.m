function [K_ic, K_ts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, alpha_ic, alpha_ts] = BankThermalProperties(Lambda)
%% calculate thermal conductivity of frozen sediment-ice mixture
% input:
%   Lambda = volume of pore space filled with ice (unitless, 0-1)
% outputs:
%   K_ic = thermal conductivity of ice cement (W/m/K)
%   K_ts = thermal conductivity of saturated, thawed bank (W/m/K)
%   rho_ic = bulk density of ice cement (kg/m3)
%   rho_ts = bulk density of mthawed sediment (kg/m3)
%   cp_ic = ice cement heat capacity (J/kg/K)
%   cp_ts = thawed bank heat capacity (J/kg/K)
%   Lf = latent heat of fusion of ice cement (J/kg)
%   alpha_ic = ice cement thermal diffusivity (m2/s)
%   alpha_ts = thawed bank thermal diffusivity (m2/s)

% constants
K_w = 0.591;                        % thermal conductivity of water (W/m/K)
cp_water = 4184;                    % specific heat of water (J/kg/K)
rho = 1000;                         % density of water (kg/m3)
K_ice = 2.18;                       % thermal conductivity of ice (W/m/K)
cp_ice = 2093;                      % specific heat of ice (J/kg/K)
rho_ice = 920;                      % density of ice (kg/m3)
Lf_ice = 333550;                    % latent heat of fusion (J/kg)
K_s = 3;                            % thermal conductivity of sediment (W/m/K)
cp_sed = 730;                       % specific heat of sediment (J/kg/K)
rho_s = 2650;                       % density of sediment (kg/m3)

% general constants
g = 9.81;                           % gravitational acceleration (m/s2)

% thermal conductivity (W/m/K)
% from "Convection in Porous Media" (2017) Section 2.2.1
K_ic = K_s^(1 - Lambda) * K_ice^Lambda;
K_ts = K_s^(1 - Lambda) * K_w^Lambda;

% bulk density (kg/m3)
rho_ic = rho_s*(1 - Lambda) + rho_ice*Lambda;
rho_ts = rho_s*(1 - Lambda) + rho*Lambda;

% % mass fraction ice and water (kg/kg)
% Mice = rho_ice*Lambda/rho_ic;
% Mw = rho*Lambda/rho_ts;

% calculate heat capacity (J/kg/K)
cp_ic = cp_ice*Lambda + cp_sed*(1-Lambda);
cp_ts = cp_water*Lambda + cp_sed*(1-Lambda);

% latent heat of fusion (J/kg)
Lf = Lf_ice * Lambda;

% thermal diffusivity (m2/s)
alpha_ic = K_ic / cp_ic / rho_ic;
alpha_ts = K_ts / cp_ts / rho_ts;


end