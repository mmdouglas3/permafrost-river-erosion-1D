% figure 9 - run model for Yukon River annual hydrograph with 10 cm bank failures

clear

% Yukon River base case

% river hydraulics
S = 1.6e-4;             % channel slope (-)
Lambda = 0.22;          % porosity (-)
Cf = 0.011;             % coefficient of friction (-)    
D_bank = 1e-3;          % bank grain size (m)
D_bed = 0.01;           % grain size (m)

% thermal properties
Tf = 0;              % temperature of fusion of water (degC)
Tbank0 = -1;         % bank background temperature (degC)

% threshold thaw layer thickness for failure (m)
yFail = 0.1;        

% unsteady model parameters
totalt = 200*60*60*24;      % model duration (s)
dt = 10;                    % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.01;                  % spatial step (m)
bankdepth = 50;             % size of bank grid (m)

% set up Yukon timeseries data
load('YukonModelParameters.mat','doy','QwFit','TwFit','Hfit','Ufit');
doyY = doy;
doy = 110 + (0:dt:totalt)/(60*60*24);
Qw = interp1(doyY,QwFit,doy);
Tw = interp1(doyY,TwFit,doy);
Tw(Tw<0) = 0;
H = Hfit.a*Qw.^Hfit.b;
U = Ufit.a*Qw.^Ufit.b;

save('YukonHydrographBF.mat');

% run unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_QwTwBF(Lambda, Cf, U', S, Tw', ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('YukonHydrograph_BankFail.mat');

