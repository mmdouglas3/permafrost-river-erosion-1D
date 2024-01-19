% run model to steady-state for ranges of parameters from global
% compilation to make figures 5 and 6

clear

% Yukon River base case

% river hydraulics
U = 1;                  % mean flow velocity (m/s)
S = 1.6e-4;             % channel slope (-)
Lambda = 0.22;          % porosity (-)
Cf = 0.011;             % coefficient of friction (-)    
D_bank = 1e-3;          % bank grain size (m)
D_bed = 0.01;           % grain size (m)
R = 1.65;               % submerged specific gravity of sediment (-)

% thermal properties
Twater = 10;         % water temperature (degC)
Tf = 0;              % temperature of fusion of water (degC)
Tbank0 = -1;         % bank background temperature (degC)

% unsteady model parameters
totalt = 60*60*72;          % model duration (s)
dt = 0.1;                   % timestep (s)
ts = totalt ./ dt;          % number of timesteps
dx = 0.01;                  % spatial step (m)
bankdepth = 20;             % size of bank grid (m)
epsilon = 1e-3;

save('YukonBaseCaseSS.mat');

%% vary lambda from base case

clear
load('YukonBaseCaseSS.mat')
disp('Vary \lambda_p (1/5)')
Lp = 0:0.05:1;

% solve unsteady model
eta_ss = zeros(size(Lp));
E = zeros(size(Lp));
M = zeros(size(Lp));
M2 = zeros(size(Lp));
eta_ts = zeros(size(Lp));
q_star = zeros(size(Lp));
for i = 1:length(Lp)-1
    [~,~, rho_ic,~,cp_ic,~, Lf,~,~] = BankThermalProperties(Lp(i));
    [TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
        RunPermafrostBankModel_TLss(Lp(i), Cf, U, S, Twater, ...
        D_bank, ts, dx, dt, bankdepth, Tbank0, epsilon);
    ind = find(round(tsP,10)<bankdepth,1,'last');
    E(i) = EP(ind);
    M(i) = MP(1);
    M2(i) = qwP(2)/(rho_ic*(Lf-Tbank0*cp_ic));
    if E(i)>0 && ~isinf(E(i))
        eta_ss(i) = mean(tsP(ind-round(100*E(i)/(dx/dt)):ind) - icP(ind-round(100*E(i)/(dx/dt)):ind));
    end
    eta_ts(i) = ind*dt;
    q_star(i) = qwP(2)/(rho_ic*Lf)/sqrt(R*9.81*D_bank);
end

% save results
save('Yukon_TLss_lambda.mat');

%% vary T_bank from base case

clear
load('YukonBaseCaseSS.mat')
disp('Vary Tbank (2/5)')
Tbank0 = -15:0.25:0;

% solve unsteady model
eta_ss = zeros(size(Tbank0));
E = zeros(size(Tbank0));
M = zeros(size(Tbank0));
M2 = zeros(size(Tbank0));
eta_ts = zeros(size(Tbank0));
q_star = zeros(size(Tbank0));
for i = 1:length(Tbank0)
    [~,~, rho_ic,~,cp_ic,~, Lf,~,~] = BankThermalProperties(Lambda);
    [TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
        RunPermafrostBankModel_TLss(Lambda, Cf, U, S, Twater, ...
        D_bank, ts, dx, dt, bankdepth, Tbank0(i), epsilon);
    ind = find(round(tsP,10)<bankdepth,1,'last');
    if isempty(ind)==true
        ind = length(tsP);
    end
    E(i) = EP(ind);
    M(i) = MP(1);
    M2(i) = qwP(2)/(rho_ic*(Lf-Tbank0(i)*cp_ic));
    if E(i)>0 && ~isinf(E(i))
        eta_ss(i) = mean(tsP(ind-round(100*E(i)/(dx/dt)):ind) - icP(ind-round(100*E(i)/(dx/dt)):ind));
    end
    eta_ts(i) = ind*dt;
    q_star(i) = qwP(2)/(rho_ic*Lf)/sqrt(R*9.81*D_bank);
end

% save results
save('Yukon_TLss_Tbank.mat');

%% vary T_water from base case

clear
load('YukonBaseCaseSS.mat')
disp('Vary Twater (3/5)')
Twater = 0:1:20;

% solve unsteady model
eta_ss = zeros(size(Twater));
E = zeros(size(Twater));
M = zeros(size(Twater));
M2 = zeros(size(Twater));
eta_ts = zeros(size(Twater));
q_star = zeros(size(Twater));
for i = 1:length(Twater)
    [~,~, rho_ic,~,cp_ic,~, Lf,~,~] = BankThermalProperties(Lambda);
    [TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
        RunPermafrostBankModel_TLss(Lambda, Cf, U, S, Twater(i), ...
        D_bank, ts, dx, dt, bankdepth, Tbank0, epsilon);
    ind = find(round(tsP,10)<bankdepth,1,'last');
    if isempty(ind)==true
        ind = length(tsP);
    end
    E(i) = EP(ind);
    M(i) = MP(1);
    M2(i) = qwP(2)/(rho_ic*(Lf-Tbank0*cp_ic));
    if E(i)>0 && ~isinf(E(i))
        eta_ss(i) = mean(tsP(ind-round(100*E(i)/(dx/dt)):ind) - icP(ind-round(100*E(i)/(dx/dt)):ind));
    end
    eta_ts(i) = ind*dt;
    q_star(i) = qwP(2)/(rho_ic*Lf)/sqrt(R*9.81*D_bank);
end

% save results
save('Yukon_TLss_Twater.mat');

%% vary U from base case

clear
load('YukonBaseCaseSS.mat')
disp('Vary U (4/5)')
U = 0.4:0.1:3;

% solve unsteady model
eta_ss = zeros(size(U));
E = zeros(size(U));
M = zeros(size(U));
M2 = zeros(size(U));
eta_ts = zeros(size(U));
q_star = zeros(size(U));
for i = 1:length(U)
    [~,~, rho_ic,~,cp_ic,~, Lf,~,~] = BankThermalProperties(Lambda);
    [TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
        RunPermafrostBankModel_TLss(Lambda, Cf, U(i), S, Twater, ...
        D_bank, ts, dx, dt, bankdepth, Tbank0, epsilon);
    ind = find(round(tsP,10)<bankdepth,1,'last');
    if isempty(ind)==true
        ind = length(tsP);
    end
    E(i) = EP(ind);
    M(i) = MP(1);
    M2(i) = qwP(2)/(rho_ic*(Lf-Tbank0*cp_ic));
    if E(i)>0 && ~isinf(E(i))
        eta_ss(i) = mean(tsP(ind-round(100*E(i)/(dx/dt)):ind) - icP(ind-round(100*E(i)/(dx/dt)):ind));
    end
    eta_ts(i) = ind*dt;
    q_star(i) = qwP(2)/(rho_ic*Lf)/sqrt(R*9.81*D_bank);
end

% save results
save('Yukon_TLss_U.mat');

%% vary D_bank from base case
clear
load('YukonBaseCaseSS.mat')
disp('Vary Dbank (5/5)')
Dbank = logspace(-5,-1,100);

% solve unsteady model
eta_ss = zeros(size(Dbank));
E = zeros(size(Dbank));
M = zeros(size(Dbank));
M2 = zeros(size(Dbank));
eta_ts = zeros(size(Dbank));
q_star = zeros(size(Dbank));
for i = 1:length(Dbank)
    [~,~, rho_ic,~,cp_ic,~, Lf,~,~] = BankThermalProperties(Lambda);
    [TP, mf, icP, tsP, qwP, MP, EP, VNstability, CFLstability] = ...
        RunPermafrostBankModel_TLss(Lambda, Cf, U, S, Twater, ...
        Dbank(i), ts, dx, dt, bankdepth, Tbank0, epsilon);
    ind = find(round(tsP,10)<bankdepth,1,'last');
    if isempty(ind)==true
        ind = length(tsP);
    end
    E(i) = EP(ind);
    M(i) = MP(1);
    M2(i) = qwP(2)/(rho_ic*(Lf-Tbank0*cp_ic));
    if E(i)>0 && ~isinf(E(i))
        eta_ss(i) = mean(tsP(ind-round(100*E(i)/(dx/dt)):ind) - icP(ind-round(100*E(i)/(dx/dt)):ind));
    end
    eta_ts(i) = ind*dt;
    q_star(i) = qwP(2)/(rho_ic*Lf)/sqrt(R*9.81*Dbank(i));
end

% save results
save('Yukon_TLss_Dbank.mat');
