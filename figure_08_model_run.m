% figure 8 - run model for Yukon River hydrograph with summer flood

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

% unsteady model parameters
totalt = 200*60*60*24;      % model duration (s)
dt = 0.1;                   % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.01;                  % spatial step (m)
bankdepth = 50;             % size of bank grid (m)

% set up Yukon timeseries data
load('YukonModelParameters.mat','doy','QwFit','TwFit','Hfit','Ufit');
doyY = doy;
doy = 110 + (0:dt:totalt)/(60*60*24);
Tw = interp1(doyY,TwFit,doy);
Tw(Tw<0) = 0;
QwRef = interp1(doyY,QwFit,doy);
Qw = fliplr(interp1(doyY,QwFit,doy)).^2 + interp1(doyY,QwFit,doy).^2;
Qw = Qw*sum(QwRef(:))/sum(Qw(:));
H = Hfit.a*Qw.^Hfit.b;
U = Ufit.a*Qw.^Ufit.b;

save('YukonHydrographSF.mat');

% run unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_QwTw(Lambda, Cf, U', S, Tw', ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save and plot results
save('YukonHydrograph_modernSF.mat');

