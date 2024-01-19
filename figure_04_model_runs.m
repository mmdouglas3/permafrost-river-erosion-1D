% permafrost unsteady bank erosion model 
% Madison Douglas (madisond@mit.edu)
% December 7, 2023
% run model to make figure 4

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
totalt = 60*60*48;          % model duration (s)
dt = 0.1;                   % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.001;                 % spatial step (m)
bankdepth = 5;              % size of bank grid (m)
% yFail = 0.001;            % threshold for thawed sediment failure (m)

save('YukonBaseCase.mat');

%% run unsteady erosion model for Yukon base case (entrainment-limited)

disp('Starting \lambda_p = 0.22 (1/7)')

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_sand_lambda22.mat');

%% run unsteady erosion model for Yukon, lambda_p = 0.5

clear
disp('Starting \lambda_p = 0.50 (2/7)')
load('YukonBaseCase.mat')
Lambda = 0.5;           % porosity (-)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_sand_lambda50.mat');

%% run unsteady erosion model for Yukon, lambda_p = 0.8

clear
disp('Starting \lambda_p = 0.80 (3/7)')
load('YukonBaseCase.mat')
Lambda = 0.8;           % porosity (-)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_sand_lambda80.mat');

%% run unsteady erosion model for Yukon, silt, lambda_p = 0.22

clear
disp('Starting silt case (4/7)')
load('YukonBaseCase.mat')
D_bank = 5e-5;              % bank grain size (m)
totalt = 60*60*72;          % model duration (s)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_CW(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_silt_lambda22.mat');

%% run unsteady erosion model for Yukon, gravel, lambda_p = 0.22

clear
disp('Starting gravel case (5/7)')
load('YukonBaseCase.mat')
D_bank = 1e-2;          % bank grain size (m)
dx = 0.01;                 % spatial step (m)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_gravel_lambda22.mat');

%% run unsteady erosion model for Yukon, gravel, lambda_p = 0.5

clear
disp('Starting \lambda_p = 0.50 (6/7)')
load('YukonBaseCase.mat')
D_bank = 1e-2;          % bank grain size (m)
Lambda = 0.5;           % porosity (-)
dx = 0.01;              % spatial step (m)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_gravel_lambda50.mat');

%% run unsteady erosion model for Yukon, gravel, lambda_p = 0.8

clear
disp('Starting \lambda_p = 0.80 (7/7)')
load('YukonBaseCase.mat')
D_bank = 1e-2;          % bank grain size (m)
Lambda = 0.8;           % porosity (-)
dx = 0.01;              % spatial step (m)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0);

% save results
save('Yukon_gravel_lambda80.mat');
