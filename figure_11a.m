% run base case to steady state to produce figure 11a

clear
load('YukonBaseCase.mat');
totalt = 210*60*60;
timesteps = dt:dt:totalt;
ts = length(timesteps);

% get bank thermal properties
[K_ic, K_ts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, alpha_ic, alpha_ts] = BankThermalProperties(Lambda);

% solve steady-state analytical equations
[Eent, Ethaw, ~, etaSS] = SteadyStateSoln(Cf, U, S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

% unsteady solution
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, ...
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

save('TL2Validation2_case1.mat');


%% run base case with U = 1.1 m/s

clear
load('YukonBaseCase.mat');
totalt = 210*60*60;
timesteps = dt:dt:totalt;
ts = length(timesteps);
U = 1.1;

% get bank thermal properties
[K_ic, K_ts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, alpha_ic, alpha_ts] = BankThermalProperties(Lambda);

% solve steady-state analytical equations
[Eent, Ethaw, ~, etaSS] = SteadyStateSoln(Cf, U, S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

% unsteady solution
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, ...
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

save('TL2Validation2_case2.mat');


%% run base case

clear
load('YukonBaseCase.mat');
totalt = 210*60*60;
timesteps = dt:dt:totalt;
ts = length(timesteps);
U = 1.2;

% get bank thermal properties
[K_ic, K_ts, rho_ic, rho_ts, cp_ic, cp_ts, Lf, alpha_ic, alpha_ts] = BankThermalProperties(Lambda);

% solve steady-state analytical equations
[Eent, Ethaw, ~, etaSS] = SteadyStateSoln(Cf, U, S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

% unsteady solution
[TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
    RunPermafrostBankModel(Lambda, Cf, U, S, ...
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

save('TL2Validation2_case3.mat');

%% plot results

clear

load('TL2Validation2_case1.mat','icP','tsP','timesteps','etaSS','totalt','D_bank');
etaSS1 = etaSS;
eta1 = movmedian(tsP-icP,1e5);
eta1 = eta1(1:1e3:end);
timesteps = timesteps(1:1e3:end);

load('TL2Validation2_case2.mat','icP','tsP','etaSS');
etaSS2 = etaSS;
eta2 = movmedian(tsP-icP,1e5);
eta2 = eta2(1:1e3:end);

load('TL2Validation2_case3.mat','icP','tsP','etaSS');
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

