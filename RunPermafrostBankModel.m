function [Tres, mfr, ic, ts, qw, M, E, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Tw, D, timesteps, dx, dt, ...
    bankdepth, Tbank0)
%% run permafrost bank erosion model
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
% outputs:
%   Tres = 11 timesteps of bank temperature through time
%   mf = melt fraction for each cell through time
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

for t = 2:timesteps
    [T, H, mf, ic(t), ts(t), qw(t), M(t), E(t)] = RunTimestep(dx, ...
        dt, ic(t-1), ts(t-1), mf, Cf, U, S, Tw, D, T, H, Lambda, Kic, Kts, rho_ic, rho_ts, cp_ic, cp_ts, Lf);
    if t == t_res(tr)
        Tres(:,tr) = T;
        mfr(:,tr) = mf;
        tr = tr + 1;
    end
    if ic(t) <= 2*dx
        disp('Timestep number: ');disp(t);
        error('Program terminated because entire ice block melted');
    end
end

% set to standard bank surface
ic = ic - dx;
ts = ts - dx;

% run stability analysis
VNstability = vonNeumann(dt, dx, [alpha_ic; alpha_ts]);
CFLstability = CourantFriedrichsLewy(max(min([M';E'])), dt, dx);

% calculate relvant timescales
Tau_E = dx/mean(E(:));
disp(['Entrainment timescale = ',num2str(Tau_E),' s']);
Tau_M = dx^2*Lf*rho_ic/(Kts*(Tw-Tbank0));
disp(['Thaw timescale = ',num2str(Tau_M),' s']);

end