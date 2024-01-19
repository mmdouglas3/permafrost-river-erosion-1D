% model validation for synthetic conditions to make figure 11b

clear

% Yukon River base case

% river hydraulics
U = 1;                  % mean flow velocity (m/s)
S = 1.6e-4;             % channel slope (-)
Lambda = 0.5;           % porosity (-)
Cf = 0.011;             % coefficient of friction (-)    
D_bank = 1e-3;          % bank grain size (m)
D_bed = 0.01;           % grain size (m)

% thermal properties
Twater = 0.1;        % water temperature (degC)
Tf = 0;             % temperature of fusion of water (degC)
Tbank0 = 0;         % bank background temperature (degC)

% unsteady model parameters
totalt = 60*60*210;         % model duration (s)
dt = 0.1;                   % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.01;                  % spatial step (m)
bankdepth = 4;              % size of bank grid (m)

% get bank thermal properties
[K_ic, K_ts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, alpha_ic, alpha_ts] = BankThermalProperties(Lambda);
K_ic = 1e2*K_ic;
K_ts = 1e2*K_ts;
alpha_ic = 1e2*alpha_ic;
alpha_ts = 1e2*alpha_ts;

save('YukonSSValidation.mat');

%% run base case

clear
load('YukonSSValidation.mat');
timesteps = dt:dt:totalt;

% solve steady-state analytical equations
[Eent, Ethaw, ~, etaSS] = SteadyStateSoln11b(Cf, U, S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

% unsteady solution
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel11b(Lambda, Cf, U, S, ...
    Twater, D_bank, ts, dx, dt, bankdepth, Tbank0);

% plot results
figure()
hold on
plot(timesteps',tsP-icP,'k-','LineWidth',2)
plot([0,totalt],etaSS*ones(1,2),'r--','LineWidth',2)
ylabel('\eta (m)')
xlabel('Time (s)')
set(gca,'FontSize',16)
box on
hold off

save('TL2Validation_case1.mat');


%% run base case with U = 1.1 m/s

clear
load('YukonSSValidation.mat');
U = 1.1;
timesteps = dt:dt:totalt;

% solve steady-state analytical equations
[Eent, Ethaw, ~, etaSS] = SteadyStateSoln11b(Cf, U, S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

% unsteady solution
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel11b(Lambda, Cf, U, S, ...
    Twater, D_bank, ts, dx, dt, bankdepth, Tbank0);

% plot results
figure()
hold on
plot(timesteps',tsP-icP,'k-','LineWidth',2)
plot([0,totalt],etaSS*ones(1,2),'r--','LineWidth',2)
ylabel('\eta (m)')
xlabel('Time (s)')
set(gca,'FontSize',16)
box on
hold off

save('TL2Validation_case2.mat');


%% run base case with U = 1.2 m/s

clear
load('YukonSSValidation.mat');
timesteps = dt:dt:totalt;
U = 1.2;

% solve steady-state analytical equations
[Eent, Ethaw, ~, etaSS] = SteadyStateSoln11b(Cf, U, S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

% unsteady solution
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel11b(Lambda, Cf, U, S, ...
    Twater, D_bank, ts, dx, dt, bankdepth, Tbank0);

% plot results
figure()
hold on
plot(timesteps',tsP-icP,'k-','LineWidth',2)
plot([0,totalt],etaSS*ones(1,2),'r--','LineWidth',2)
ylabel('\eta (m)')
xlabel('Time (s)')
set(gca,'FontSize',16)
box on
hold off

save('TL2Validation_case3.mat');

%% plot results

clear

load('TL2Validation_case1.mat','icP','tsP','timesteps','etaSS','totalt','D_bank');
etaSS1 = etaSS;
eta1 = movmedian(tsP-icP,1e5);
eta1 = eta1(1:1e3:end);
timesteps = timesteps(1:1e3:end);

load('TL2Validation_case2.mat','icP','tsP','etaSS');
etaSS2 = etaSS;
eta2 = movmedian(tsP-icP,1e5);
eta2 = eta2(1:1e3:end);

load('TL2Validation_case3.mat','icP','tsP','etaSS');
etaSS3 = etaSS;
eta3 = movmedian(tsP-icP,1e5);
eta3 = eta3(1:1e3:end);

figure()
hold on
plot(timesteps'/60/60,eta1/D_bank,'-','LineWidth',2)
plot([0,totalt]/60/60,etaSS1*ones(1,2)/D_bank,'--','LineWidth',3)
plot(timesteps'/60/60,eta2/D_bank,'-','LineWidth',2)
plot([0,totalt]/60/60,etaSS2*ones(1,2)/D_bank,'--','LineWidth',3)
plot(timesteps'/60/60,eta3/D_bank,'-','LineWidth',2)
plot([0,totalt]/60/60,etaSS3*ones(1,2)/D_bank,'--','LineWidth',3)
legend({'numerical, {\itU} = 1.0 m/s','analytical, {\itU} = 1.0 m/s', ...
    'numerical, {\itU} = 1.1 m/s','analytical, {\itU} = 1.1 m/s', ...
    'numerical, {\itU} = 1.2 m/s','analytical, {\itU} = 1.2 m/s'},'location','southeast')
xlim([0,200])
ylabel('{\it\eta*}')
xlabel('Time (hours)')
set(gca,'FontSize',16)
box on
hold off

