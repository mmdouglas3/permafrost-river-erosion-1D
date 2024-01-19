function [Eent, Ethaw, qw, etaSS, Eent_star, Tw_star, U_star, T0_star, Cic_star, Ch_star, Kts2_star] = ...
    SteadyStateSoln(Cf, U, S, D, Twater, Lambda, Tbank0, Tf, s)

% bank thermal properties and constants
rho = 1000;         % density of water (kg/m3)
R = 1.65;           % submerged specific gravity of water (-)
g = 9.81;           % gravitational acceleration (m/s2)
cp_w = 4184;        % specific heat of water (J/kg/K)
[~, Kts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, ~, ~] = BankThermalProperties(Lambda);

% calculate bank erosion rate (m/s)
H = Cf * U^2 / (9.81 * S);                      % normal flow depth (m)
ks_b = 2.5*(2.2*D);                             % bank roughness coefficient (m)
[Cf_b,ustar] = ManningStrickler(U, S, ks_b);    % bank coefficient of friction (-)
if ks_b*ustar*1e6 < 100
    [~,Cf_b] =  CWDarcyWeisbach(H, U, ks_b);
end
Tau_crit = Parker2003Shields(D)*1650*g*D;       % critical shear stress for entrainment (Pa)
Eent = Partheneides1965(U, Cf_b, Tau_crit, 2650, 1-Lambda, 0.25e-4);

% calculate bank thaw rate (m/s)
[qw, Ch] = YaglomKader1974_transitional(U, H, ks_b, Cf_b, Twater, Tbank0);
Ethaw = qw / (rho_ic * (Lf + cp_ic*(Tf - Tbank0)));

% calculate nondimensional parameters           
Kts2_star = Kts / (rho_ts*cp_ts*sqrt(R*g*D^3));
Ch_star = (rho*cp_w*Ch) / (rho_ic*cp_ic);
U_star = U / sqrt(R*g*D);
Tw_star = (Twater*cp_ic) / Lf;
Eent_star = Eent / sqrt(R*g*D);
T0_star = (Tbank0*cp_ic) / Lf;
Cic_star = (rho_ic*cp_ic) / (rho_ts*cp_ts);

% calculate thaw layer thickness (m)
etaSS =  D*(s*Kts2_star*(Ch_star*U_star*Tw_star - Eent_star*(1-T0_star))) / ...
    (Eent_star*(1-T0_star)*(Cic_star*Ch_star*U_star + Eent_star));

end