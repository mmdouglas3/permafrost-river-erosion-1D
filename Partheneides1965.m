function E = Partheneides1965(U, Cf_b, Tau_crit, rho_sed, f_sed, M)
%% calculate entrainment-limited riverbank erosion for cohesive sediment under normal flow conditions using the entrainment equation from Partheneides (1965)
% inputs:
%   U = water velocity (m/s)
%   Cf_b = bank coefficient of friction (-)
%   Tau_crit = critical shear stress to entrain sediment (Pa)
%   rho_sed = bulk density of riverbank sediment (kg/m3)
%   f_sed = volumetric fraction of riverbank comprised of sediment (m3/m3)
%   M = coefficient for entrainment equation (kg/m2/s)
% outputs:
%   E = entrainment-limited bank erosion rate (m/s)

% constants
rho_w = 1000;       % water density (kg/m3)
n = 1;              % exponent for entrainment equation (dimensionless)

% calculate shear stress on the riverbank
Tau_bank = rho_w*Cf_b*U^2;

% if bank stress is greater than the threshold for entrainment, calculate
% erosion rates
if Tau_bank > Tau_crit
    E = M*(Tau_bank / Tau_crit - 1)^n;
else
    E = 0;
end

% convert erosion rate from kg/m2/s to m/s
E = E / (rho_sed*f_sed);

end