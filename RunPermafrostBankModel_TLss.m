function [Tres, mfr, ic, ts, qw, M, E, VNstability, CFLstability] = ...
    RunPermafrostBankModel_TLss(Lambda, Cf, U, S, Tw, D, timesteps, dx, ...
    dt, bankdepth, Tbank0, epsilon)
%% run permafrost bank erosion model to steady-state thaw layer thickness
% inputs:
%   Lambda = bed porosity (0 - 1)
%   Cf = open channel coefficient of friction (-)
%   U = mean open channel flow velocity (m/s)
%   S = channel slope (m/m)
%   Tw = average water temperature (degC)
%   D = median grain size (m)
%   timesteps = number of timesteps to run model (-)
%   dx = grid spacing (m)
%   dt = timestep spacing (s)
%   bankdepth = bankfull depth of open channel (m)
%   Tbank0 = initial temperature of frozen bank (degC)
%   epsilon = fraction change in thaw layer thickness for steady-state
% outputs:
%   Tres = 11 timesteps of bank temperature through time
%   mfr = melt fraction for each cell through time
%   ic = location of ice cement thaw front through time
%   ts = location of erosion front through time
%   qw = heat flux to wall through time (J/m2)
%   M = melt rate through time (m/s)
%   E = erosion rate through time (m/s)
%   VNstability = t/f if passes von Neumann stability criterion
%   CFLstability = t/f if passes CFL stability criterion

[Kic, Kts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, alpha_ic, alpha_ts] = BankThermalProperties(Lambda);

x = transpose(0:dx:bankdepth);                  % distance from wall (m)
mf = zeros(length(x),1);                        % fraction of ice thawed per cell (-)
T = Tbank0*ones(length(x)+1,1);                 % temperature (degC)
T(end,:) = Tw;                                  % water boundary layer (degC)
H = rho_ic*cp_ic*T*dx;                          % enthalpy (J)
H(end) = (rho_ts*cp_ts*Tw + rho_ic*Lf)*dx;             % enthalpy for water boundary layer (J)
ic = (bankdepth+dx)*ones(int64(timesteps),1);          % ice cement thaw front (m)
ts = (bankdepth+dx)*ones(int64(timesteps),1);          % thawed sediment erosion front (m)
qw = zeros(int64(timesteps),1);                        % heat flux to bank (J/m2/s)
M = zeros(int64(timesteps),1);                         % thaw rate of ice cement (m/s)
E = zeros(int64(timesteps),1);                         % erosion rate of bank (m/s)

mfr = zeros(length(mf),10);                     % results of melt fraction
mfr(:,1) = mf;                                  % initial melt fraction
Tres = zeros(length(T),10);                     % results of timesteps
Tres(:,1) = T;                                  % initial timesteps
tr = 1;                                         % save counter
t_res = transpose(timesteps/10:timesteps/10:timesteps);     % times to save T profile
t = 2;
dL = 1;

% calculate bank coefficient of friction
ks_b = 2.5*(2.2*D);                             % bank roughness coefficient (m)
Cf_b = ManningStrickler(U, S, ks_b);            % bank coefficient of friction (-)
Re_ks = U*sqrt(Cf_b) * ks_b * 1e6;              % bank roughness Reynolds number (-)
if Re_ks < 200
    Cf_b1 = Cf_b;
    [~, Cf_b2] = CWDarcyWeisbach(Cf*U^2/(9.81*S), U, ks_b);
    Cf_b = max([Cf_b1;Cf_b2]);
end

while abs(dL) > epsilon && t < timesteps
    [T, H, mf, ic(t), ts(t), qw(t), M(t), E(t)] = RunTimestep_TLss(dx, ...
        dt, ic(t-1), ts(t-1), mf, Cf, U, S, Tw, D, T, H, Lambda, Kic, Kts, ...
        rho_ic, rho_ts, cp_ic, cp_ts, Lf, Cf_b);
    % save results at interval tr
    if t == t_res(tr)
        Tres(:,tr) = T;
        mfr(:,tr) = mf;
        tr = tr + 1;
    end
    % error message if bank too small
    if ic(t) <= 2*dx
        disp('Timestep number: ');disp(t);
        error('Program terminated because entire ice block melted');
    end
    % determine if steady-state reached
    if ic(t) < bankdepth-2*dx && t > round(100*E(t)/(dx/dt))
        dL = (mean(ts(t-round(100*E(t)/(dx/dt)):t) - ic(t-round(100*E(t)/(dx/dt)):t)) - ts(t) + ic(t)) ...
            / mean(ts(t-round(100*E(t)/(dx/dt)):t) - ic(t-round(100*E(t)/(dx/dt)):t));
    end
    t = t + 1;
end

% set to standard bank surface
ic = ic - dx;
ts = ts - dx;

% run stability analysis
VNstability = vonNeumann(dt, dx, [alpha_ic; alpha_ts]);
CFLstability = CourantFriedrichsLewy(max(min([M';E'])), dt, dx);

end