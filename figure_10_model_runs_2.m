% run model for 1 cm sediment to make figure 10

clear

% Yukon River base case

% river hydraulics
U = 1;                  % mean flow velocity (m/s)
S = 1.6e-4;             % channel slope (-)
Lambda = 0.22;          % porosity (-)
Cf = 0.011;             % coefficient of friction (-)    
D_bank = 1e-2;          % bank grain size (m)
D_bed = 0.01;           % grain size (m)

% thermal properties
Twater = 15;         % water temperature (degC)
Tf = 0;              % temperature of fusion of water (degC)
Tbank0 = -1;         % bank background temperature (degC)

% unsteady model parameters
totalt = 60*60*24;          % model duration (s)
dt = 0.1;                   % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.01;                  % spatial step (m)
bankdepth = 15;             % size of bank grid (m)

save('YukonFailCase.mat');

%% run unsteady erosion model for Yukon, 0.1 cm bank failure threshold

disp('Starting failure threshold = 0.1 cm (1/8)')
yFail = 0.001;          % bank threshold for failure (m)

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_o1cmFail_Tw15C.mat');

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
save('Yukon_nofail_Tw15C.mat');

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
save('Yukon_1cmFail_Tw15C.mat');

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
save('Yukon_3cmFail_Tw15C.mat');

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
save('Yukon_5cmFail_Tw15C.mat');

%% run unsteady erosion model for Yukon, 10 cm bank failure

clear
disp('Starting failure threshold = 10 cm (6/8)')
load('YukonFailCase.mat')
yFail = 0.10;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_10cmFail_Tw15C.mat');

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
save('Yukon_15cmFail_Tw15C.mat');

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
save('Yukon_25cmFail_Tw15C.mat');

%% run unsteady erosion model for Yukon, 0.1 cm bank failure threshold

clear
disp('Starting failure threshold = 0.1 cm (1/8)')
load('YukonFailCase.mat')
yFail = 0.001;          % bank threshold for failure (m)
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_o1cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, no bank failure

clear
disp('No bank failure (2/8)')
load('YukonFailCase.mat')
yFail = 30;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_nofail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 1 cm bank failure

clear
disp('Starting failure threshold = 1 cm (3/8)')
load('YukonFailCase.mat')
yFail = 0.01;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_1cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 3 cm bank failure

clear
disp('Starting failure threshold = 3 cm (4/8)')
load('YukonFailCase.mat')
yFail = 0.03;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_3cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 5 cm bank failure

clear
disp('Starting failure threshold = 5 cm (5/8)')
load('YukonFailCase.mat')
yFail = 0.5;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_5cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 10 cm bank failure

clear
disp('Starting failure threshold = 10 cm (6/8)')
load('YukonFailCase.mat')
yFail = 0.10;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_10cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 15 cm bank failure

clear
disp('Starting failure threshold = 15 cm (7/8)')
load('YukonFailCase.mat')
yFail = 0.15;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_15cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 25 cm bank failure

clear
disp('Starting failure threshold = 25 cm (8/8)')
load('YukonFailCase.mat')
yFail = 0.25;
Twater = 5;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_25cmFail_Tw5C.mat');

%% run unsteady erosion model for Yukon, 0.1 cm bank failure threshold

clear
disp('Starting failure threshold = 0.1 cm (1/8)')
load('YukonFailCase.mat')
yFail = 0.001;          % bank threshold for failure (m)
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_o1cmFail_Tw1C.mat');

%% run unsteady erosion model for Yukon, no bank failure

clear
disp('No bank failure (2/8)')
load('YukonFailCase.mat')
yFail = 30;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_nofail_Tw1C.mat');

%% run unsteady erosion model for Yukon, 1 cm bank failure

clear
disp('Starting failure threshold = 1 cm (3/8)')
load('YukonFailCase.mat')
yFail = 0.01;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_1cmFail_Tw1C.mat');

%% run unsteady erosion model for Yukon, 3 cm bank failure

clear
disp('Starting failure threshold = 3 cm (4/8)')
load('YukonFailCase.mat')
yFail = 0.03;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_3cmFail_Tw1C.mat');

%% run unsteady erosion model for Yukon, 5 cm bank failure

clear
disp('Starting failure threshold = 5 cm (5/8)')
load('YukonFailCase.mat')
yFail = 0.05;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_5cmFail_Tw1C.mat');

%% run unsteady erosion model for Yukon, 10 cm bank failure

clear
disp('Starting failure threshold = 10 cm (6/8)')
load('YukonFailCase.mat')
yFail = 0.10;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_10cmFail_Tw1C.mat');

%% run unsteady erosion model for Yukon, 15 cm bank failure

clear
disp('Starting failure threshold = 15 cm (7/8)')
load('YukonFailCase.mat')
yFail = 0.15;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_15cmFail_Tw1C.mat');

%% run unsteady erosion model for Yukon, 25 cm bank failure

clear
disp('Starting failure threshold = 25 cm (8/8)')
load('YukonFailCase.mat')
yFail = 0.25;
Twater = 1;

% solve unsteady model
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel_Fail(Lambda, Cf, U, S, Twater, ...
    D_bank, ts, dx, dt, bankdepth, Tbank0, yFail);

% save results
save('Yukon_25cmFail_Tw1C.mat');

