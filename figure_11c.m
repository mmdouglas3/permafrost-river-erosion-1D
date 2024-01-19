% make plot of steady state thaw layer thickness versus the analytical
% solution fitted to temperature profile factor s

clear

%% get steady state conditions for Dbank

load('YukonBaseCaseSS.mat');
load('Dbank_SS.mat','eta_ts','Dbank');
eta_Dbank = eta_ts;

% analytical steady-state solution for each set of conditions
for i = 1:length(Dbank)
    [Eent, ~, ~, eta_Dbank_ss(i), Dbank_Eents(i), Dbank_Tws(i), Dbank_Us(i), ...
        Dbank_T0s(i), Dbank_Cics(i), Dbank_Chs(i), Dbank_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, Dbank(i), Twater, Lambda, Tbank0, Tf, 1);
    [~,~,~,~,~,~,~,~, alphats] = BankThermalProperties(Lambda);
    R_Dbank(i) = alphats ./ Eent ./ eta_Dbank_ss(i);

    [~,~,~, etaSS_Dbank_fit(i)] = SteadyStateSoln_fittedS(Cf, U, S, Dbank(i), Twater, Lambda, Tbank0, Tf);
end


%% get steady state conditions for Lambda

load('YukonBaseCaseSS.mat');
load('Lambda_SS.mat','eta_ts','Lp');
eta_Lambda = eta_ts;

% analytical solution for each set of conditions
for i = 1:length(Lp)
    [Eent, ~, ~, eta_Lambda_ss(i), Lambda_Eents(i), Lambda_Tws(i), Lambda_Us(i), ...
        Lambda_T0s(i), Lambda_Cics(i), Lambda_Chs(i), Lambda_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, D_bank, Twater, Lp(i), Tbank0, Tf, 1);
    [~,~,~,~,~,~,~,~, alphats] = BankThermalProperties(Lp(i));
    R_Lambda(i) = alphats ./ Eent ./ eta_Lambda_ss(i);

    [~,~,~, etaSS_Lambda_fit(i)] = SteadyStateSoln_fittedS(Cf, U, S, D_bank, Twater, Lp(i), Tbank0, Tf);
end

%% get steady state conditions for Tbank0

load('YukonBaseCaseSS.mat');
load('Tbank0_SS.mat','eta_ts','Tbank0');
eta_Tbank0 = eta_ts;
[~,~,~,~,~,~,~,~, alphats] = BankThermalProperties(Lambda);

% analytical solution for each set of conditions
for i = 1:length(Tbank0)
    [Eent, ~, ~, eta_Tbank0_ss(i), Tbank0_Eents(i), Tbank0_Tws(i), Tbank0_Us(i), ...
        Tbank0_T0s(i), Tbank0_Cics(i), Tbank0_Chs(i), Tbank0_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, D_bank, Twater, Lambda, Tbank0(i), Tf, 1);
    R_Tbank0(i) = alphats ./ Eent ./ eta_Tbank0_ss(i);

    [~,~,~, etaSS_Tbank0_fit(i)] = SteadyStateSoln_fittedS(Cf, U, S, D_bank, Twater, Lambda, Tbank0(i), Tf);
end


%% get steady state conditions for Twater

load('YukonBaseCaseSS.mat');
load('Twater_SS.mat','eta_ts','Twater');
eta_Twater = eta_ts;

% analytical solution for each set of conditions
for i = 1:length(Twater)
    [Eent, ~, ~, eta_Twater_ss(i), Twater_Eents(i), Twater_Tws(i), Twater_Us(i), ...
        Twater_T0s(i), Twater_Cics(i), Twater_Chs(i), Twater_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, D_bank, Twater(i), Lambda, Tbank0, Tf, 1);
    R_Twater(i) = alphats ./ Eent ./ eta_Twater_ss(i);

    [~,~,~, etaSS_Twater_fit(i)] = SteadyStateSoln_fittedS(Cf, U, S, D_bank, Twater(i), Lambda, Tbank0, Tf);
end


%% get steady state conditions for U

load('YukonBaseCaseSS.mat');
load('U_SS.mat','eta_ts','U');
eta_U = eta_ts;

% analytical solution for each set of conditions
for i = 1:length(U)
    [Eent, ~, ~, eta_U_ss(i), U_Eents(i), U_Tws(i), U_Us(i), ...
        U_T0s(i), U_Cics(i), U_Chs(i), U_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U(i), S, D_bank, Twater, Lambda, Tbank0, Tf, 1);
    R_U(i) = alphats ./ Eent ./ eta_U_ss(i);

    [~,~,~, etaSS_U_fit(i)] = SteadyStateSoln_fittedS(Cf, U(i), S, D_bank, Twater, Lambda, Tbank0, Tf);
end

%% Try fit logarithmic regression

y = [eta_Dbank./eta_Dbank_ss, eta_Lambda./eta_Lambda_ss, ...
    eta_Tbank0./eta_Tbank0_ss, eta_Twater./eta_Twater_ss, eta_U./eta_U_ss]';
y2 = log10(y);
ymask = ~isinf(y2) & ~isnan(y2) & imag(y2)==0;
ymask(182) = 0;
y = y2(ymask);

Eents = [Dbank_Eents, Lambda_Eents, Tbank0_Eents, Twater_Eents, U_Eents]';
T0s = [Dbank_T0s, Lambda_T0s, Tbank0_T0s, Twater_T0s, U_T0s]';
Us = [Dbank_Us, Lambda_Us, Tbank0_Us, Twater_Us, U_Us]';
Chs = [Dbank_Chs, Lambda_Chs, Tbank0_Chs, Twater_Chs, U_Chs]';
Cics = [Dbank_Cics, Lambda_Cics, Tbank0_Cics, Twater_Cics, U_Cics]';
Tws = [Dbank_Tws, Lambda_Tws, Tbank0_Tws, Twater_Tws, U_Tws]';
Kts2s = [Dbank_Kts2s, Lambda_Kts2s, Tbank0_Kts2s, Twater_Kts2s, U_Kts2s]';
X1 = [Eents, abs(T0s), Us, Chs, Cics, Tws, Kts2s];
X1 = X1(ymask,:);
X2 = log10(X1);
Sfunc = fitlm(X2,y);

save('FittedSfunction.mat','Sfunc');
%% plot comparison between fit and actual solution

figure()
hold on
plot([1e2,1e6],[1e2,1e6],'k-');
plot(eta_Dbank, etaSS_Dbank_fit,'o','LineWidth',2);
plot(eta_Lambda, etaSS_Lambda_fit,'ks','LineWidth',2);
plot(eta_Tbank0, etaSS_Tbank0_fit,'d','LineWidth',2);
plot(eta_Twater, etaSS_Twater_fit,'>','LineWidth',2);
plot(eta_U, etaSS_U_fit,'*','LineWidth',2);
xlabel('{\it\eta_s_s}');
ylabel('{\it\eta_{ss,a}}');
set(gca,'FontSize',16,'YScale','log','XScale','log');
legend({'1:1 line','Vary {\itD_5_0}','Vary {\it\lambda_p}','Vary {\itT_{bank}}', ...
    'Vary {\itT_{water}}','Vary {\itU}'},'location','southeast');
box on
grid on
hold off
