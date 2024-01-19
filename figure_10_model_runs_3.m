% run model for 1 mm sediment to make figure 10

clear

% Yukon River base case

% river hydraulics
U = 1;                  % mean flow velocity (m/s)
S = 1.6e-4;             % channel slope (-)
Lambda = 0.22;          % porosity (-)
Cf = 0.011;             % coefficient of friction (-)    
D_bank = 1e-3;          % bank grain size (m)
D_bed = 0.01;           % grain size (m)

% thermal properties
Twater = 10;         % water temperature (degC)
Tf = 0;              % temperature of fusion of water (degC)
Tbank0 = -1;         % bank background temperature (degC)

% unsteady model parameters
totalt = 60*60*24;          % model duration (s)
dt = 0.1;                   % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.01;                  % spatial step (m)
bankdepth = 15;             % size of bank grid (m)

save('YukonFailCase.mat');

%% run unsteady erosion model for Yukon, 1 mm bank failure threshold

disp('Starting failure threshold = 1 mm (1/8)')
yFail = 0.001;          % bank threshold for failure (m)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_1mmFail.mat');

%% run unsteady erosion model for Yukon, no bank failure

clear
disp('No bank failure (2/8)')
load('YukonFailCase.mat')
yFail = 30;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_nofail_D1mm.mat');

%% run unsteady erosion model for Yukon, 1 cm bank failure

clear
disp('Starting failure threshold = 1 cm (3/8)')
load('YukonFailCase.mat')
yFail = 0.01;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_10mmFail.mat');

%% run unsteady erosion model for Yukon, 3 cm bank failure

clear
disp('Starting failure threshold = 3 cm (4/8)')
load('YukonFailCase.mat')
yFail = 0.03;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_30mmFail.mat');

%% run unsteady erosion model for Yukon, 5 cm bank failure

clear
disp('Starting failure threshold = 5 cm (5/8)')
load('YukonFailCase.mat')
yFail = 0.05;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_50mmFail.mat');

%% run unsteady erosion model for Yukon, 10 cm bank failure

clear
disp('Starting failure threshold = 10 cm (6/8)')
load('YukonFailCase.mat')
yFail = 0.1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_100mmFail.mat');

%% run unsteady erosion model for Yukon, 15 cm bank failure

clear
disp('Starting failure threshold = 15 cm (7/8)')
load('YukonFailCase.mat')
yFail = 0.15;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_150mmFail.mat');

%% run unsteady erosion model for Yukon, 25 cm bank failure

clear
disp('Starting failure threshold = 25 cm (8/8)')
load('YukonFailCase.mat')
yFail = 0.25;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_250mmFail.mat');

