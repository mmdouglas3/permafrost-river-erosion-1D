function [Tbank2, Hbank2, mf2, ic2, ts2, qw, M, E] = RunTimestep_TLss(dx, ...
    dt, ic, ts, mf, Cf, U, S, Tw, D, Tbank, Hbank, Lambda, Kic, Kts, rho_ic, ...
    rho_ts, cp_ic, cp_ts, Lf, Cf_b)
%% run timestep of bank thaw then erosion
% inputs:
%   dx = spatial step size (m)
%   dt = timestep size (s)
%   ic = location of ice cement thaw front (m)
%   ts = location of thawed erosional front (m)
%   mf = fraction of ice melted in each cell (0 - 1)
%   Cf = coefficient of friction in channel (-)
%   U = mean flow velocity in channel (m/s)
%   S = channel slope (m/m)
%   Tw = mean temperature of channel water (degC)
%   D = median grain size (m)
%   Tbank = vector of bank temperatures through space (degC)
%   Hbank = enthalpy of bank (temp + latent heat) (J)
%   Lambda = bank porosity (vol ice / vol sed)
%   Kic = thermal conductivity of ice cement (W/m/K)
%   Kts = thermal conductivity of saturated, thawed bank (W/m/K)
%   rho_ic = bulk density of ice cement (kg/m3)
%   rho_ts = bulk density of mthawed sediment (kg/m3)
%   cp_ic = ice cement heat capacity (J/kg/K)
%   cp_ts = thawed bank heat capacity (J/kg/K)
%   Lf = latent heat of fusion of ice cement (J/kg)
%   Cf_b = bank coefficient of friction (-)
% outputs:
%   Tbank2 = updated temperature profile in bank (degC)
%   Hbank2 = updated enthalpy of bank (temp + latent heat) (J)
%   mf2 = updated fraction of ice melted in each cell (0 - 1)
%   ic2 = new location of ice cement thaw front (m)
%   ts2 = new location of thawed sediment erosion front (m)
%   qw = heat flux to bank (J/m2/s)
%   M = melt rate (m/s)
%   E = erosion rate (m/s)

% find location of melt and erosion fronts
mf_ind = ceil(round(ic,10)/dx);
ef_ind = ceil(round(ts,10)/dx);
if mf_ind > length(mf)
    mf_ind = length(mf);
    ef_ind = mf_ind;
end
ef = ef_ind - ts/dx;

% calculate bank erosion rate (m/s)
H = Cf * U^2 / (9.81 * S);                     % normal flow depth (m)
Tau_crit = Parker2003Shields(D)*1650*9.81*D;   % critical shear stress for entrainment (Pa)
E = Partheneides1965(U, Cf_b, Tau_crit, 2650, 1-Lambda, 0.25e-4);

% assign thermal properties throughout bank
Kbank = Kic*ones(size(Tbank(1:ef_ind)));

% calculate thermal properties for thawed layer
if ic < ts
    % calculate thawed layer flow velocities (m/s)                
    [K, ~, ~] = LapotreLamb2018(D);         % average sediment permeability (m2)
    [Usub1, ~, P] = Lamb2017WRR_limited(H, S, D, K, Lambda);
    
    % find location of fluid velocity boundary
    P_ind = ceil(round(ts-P,10)/dx);
    
    % calculate thermal conductivity (W/m/s) = molec + dispersion
    Kt1 = HsuCheng1990(Usub1, D);
    Kts1 = Kt1 + Kts;
    Kts2 = Kts;

    % 1) thaw depth less than exchange layer thickness, melt front and
    % erosion front in different pixels
    if ts-P < ic && mf_ind < ef_ind
        Kbank(mf_ind+1:end) = Kts1;
        Kbank(mf_ind) = (mf(mf_ind)/Kts1 + (1-mf(mf_ind))/Kic)^(-1);
        
    % 2) thaw depth less than exchange layer thickness, melt front and
    % erosion front in same pixel
    elseif ts-P < ic && mf_ind == ef_ind
        f_ts1 = mf(mf_ind)/(1-ef);
        f_ic = 1 - f_ts1;
        Kbank(mf_ind) = (f_ts1/Kts1 + f_ic/Kic)^(-1);

    % 3) exchange layer boundary within melt front pixel, erosion front in
    % different pixel
    elseif P_ind == mf_ind && mf_ind < ef_ind
        Kbank(mf_ind) = ((mf(mf_ind) - rem(ts-P,dx))/Kts1 + rem(ts-P,dx)/Kts2 + (1-mf(mf_ind))/Kic)^(-1);
        Kbank(P_ind+1:end) = Kts1;
        
    % 4) exchange layer boundary within melt front pixel, erosion front in
    % same pixel
    elseif P_ind == mf_ind && mf_ind == ef_ind
        f_ic = (1-mf(mf_ind))/(1-ef);
        f_ts1 = P/(1-ef);
        f_ts2 = 1 - f_ic - f_ts1;
        Kbank(mf_ind) = (f_ts1/Kts1 + f_ts2/Kts2 + f_ic/Kic)^(-1);

    % 5) exchange layer boundary within erosion front pixel, melt front in
    % different pixel
    elseif P_ind == ef_ind && mf_ind < ef_ind
        f_ts1 = P/(1-ef);
        f_ts2 = 1 - f_ts1;
        Kbank(P_ind) = (f_ts1/Kts1 + f_ts2/Kts2)^(-1);
        Kbank(mf_ind) = (mf(mf_ind)/Kts2 + (1-mf(mf_ind))/Kic)^(-1);
        if P_ind > mf_ind+1
            Kbank(mf_ind+1:P_ind-1) = Kts2;
        end

    % 6) melt front, erosion front, and exchange layer boundaries all in
    % different pixels
    elseif P_ind < ef_ind && P_ind > mf_ind
        Kbank(mf_ind) = (mf(mf_ind)/Kts2 + (1-mf(mf_ind))/Kic)^(-1);
        Kbank(P_ind) = (rem(ts-P, dx)/Kts2 + (1-rem(ts-P, dx))/Kts1)^(-1);
        Kbank(P_ind+1:ef_ind) = Kts1;
        if P_ind > mf_ind+1
            Kbank(mf_ind+1:P_ind-1) = Kts2;
        end
    end
end

% calculate and distribute heat flux to bank
ks_b = 2.5*(2.2*D);                             % bank roughness coefficient (m)
Tb = Tbank(ef_ind);
qw = YaglomKader1974_transitional(U, H, ks_b, Cf_b, Tw, Tb);    % heat flux to bank (J/m2/s)
Hbank2 = Hbank;                                             % bank enthalpy (J)
Hbank2(ef_ind) = Hbank2(ef_ind) + qw*dt;

if Hbank2(ef_ind) > rho_ic*Lf*dx
    Tbank(ef_ind) = (Hbank2(ef_ind) - rho_ic*Lf*dx)/(dx*rho_ts*cp_ts);
elseif Hbank2(ef_ind) > 0
    Tbank(ef_ind) = 0;
else
    Tbank(ef_ind) = Hbank2(ef_ind) / (dx*rho_ic*cp_ic);
end

% update bank enthalpy (J/m2) using upwinding finite differences
dH = zeros(size(Hbank2));
dH(1) = dt*Kbank(2)*(Tbank(2)-Tbank(1))/dx;
% no heat removed from last grid node, may cause fast thaw
dH(2:ef_ind-1) = dt*Kbank(3:ef_ind).*diff(Tbank(2:ef_ind))/dx ...
    - dt*Kbank(2:ef_ind-1).*diff(Tbank(1:ef_ind-1))/dx;
dH(ef_ind) = -dt*Kbank(ef_ind)*(Tbank(ef_ind)-Tbank(ef_ind-1))/dx;
Hbank2 = Hbank2 + dH;

% update melt fraction
mf2 = Hbank2(1:end-1)/(rho_ic*Lf*dx);
mf2(mf2<0) = 0;
mf2(mf2>1) = 1;
% calculate melt front movement rate (m/s)
M = (sum(mf2) - sum(mf))*dx/dt;

% update bank temperature
Tbank2 = zeros(size(Tbank));
Tbank2(mf2 == 0) = Hbank2(mf2 == 0)/(dx*rho_ic*cp_ic);
Tbank2(mf2 == 1) = (Hbank2(mf2 == 1) - rho_ic*Lf*dx)/(dx*rho_ts*cp_ts);

% calculate new erosional and melt fronts
ic2 = ic - M*dt;           % ice cement melt front (m)
ts2 = ts - E*dt;           % sediment erosional front (m)
if ts2 < ic2
    ts2 = ic2;          % frozen deposits can't be eroded
end

% set all eroded cells to have open channel water temperature
ef_ind2 = ceil(round(ts2,10)/dx);
Tbank2(ef_ind2+1:end) = Tw;
Hbank2(ef_ind2+1:end) = (rho_ts*cp_ts*Tw + rho_ic*Lf)*dx;

end