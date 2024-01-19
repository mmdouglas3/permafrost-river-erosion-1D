function Keff = HsuCheng1990(U, D50)
%% calculate effective conductivity from thermal dispersion in a saturated porous medium
% inputs:
%   U = mean velocity within subsurface layer (m/s)
%   D50 = length scale (m), typically median bank grain size or sqrt(hydraulic permittivity)
% output:
%   Keff = effective thermal dispersion conductivity (W/m/K)

% load water thermal properties
K_w = 0.591;                        % thermal conductivity of water (W/m/K)
alpha_w = 1.32*10^(-7);             % molecular diffusivity (m2/s)

% effective thermal conductivity due to thermal dispersion (W/m/K)
Dt = 0.12;                      % constant from paper (unitless)
Keff = K_w * Dt * (D50 / alpha_w) * U;

end