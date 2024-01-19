% make plot of steady state thaw layer thickness versus the analytical
% solution as a function of their mis-match

clear

% get steady state conditions for Dbank
load('YukonBaseCaseSS.mat');
load('Yukon_TLss_Dbank.mat','eta_ss','Dbank','eta_ts','M','E');
tss_Dbank = eta_ts;
tssa_Dbank = eta_ss./(M-E);
for i = 1:length(Dbank)
    [Eent, ~, ~, eta_Dbank_ss(i), Dbank_Eents(i), Dbank_Tws(i), Dbank_Us(i), ...
        Dbank_T0s(i), Dbank_Cics(i), Dbank_Chs(i), Dbank_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, Dbank(i), Twater, Lambda, Tbank0, Tf, 1);

    [~,~,~, tssa_Dbank_fit(i)] = SteadyStateSoln_fittedTss(Cf, U, S, Dbank(i), Twater, Lambda, Tbank0, Tf);
end

% get steady state conditions for Lambda
load('YukonBaseCaseSS.mat');
load('Yukon_TLss_lambda.mat','eta_ss','Lp','eta_ts','M','E');
tss_Lambda = eta_ts;
tssa_Lambda = eta_ss./(M-E);
for i = 1:length(Lp)
    [Eent, ~, ~, eta_Lambda_ss(i), Lambda_Eents(i), Lambda_Tws(i), Lambda_Us(i), ...
        Lambda_T0s(i), Lambda_Cics(i), Lambda_Chs(i), Lambda_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, D_bank, Twater, Lp(i), Tbank0, Tf, 1);

    [~,~,~, tssa_Lambda_fit(i)] = SteadyStateSoln_fittedTss(Cf, U, S, D_bank, Twater, Lp(i), Tbank0, Tf);
end

% get steady state conditions for Tbank0
load('YukonBaseCaseSS.mat');
load('Yukon_TLss_Tbank.mat','eta_ss','Tbank0','eta_ts','M','E');
tss_Tbank0 = eta_ts;
tssa_Tbank0 = eta_ss./(M-E);
for i = 1:length(Tbank0)
    [Eent, ~, ~, eta_Tbank0_ss(i), Tbank0_Eents(i), Tbank0_Tws(i), Tbank0_Us(i), ...
        Tbank0_T0s(i), Tbank0_Cics(i), Tbank0_Chs(i), Tbank0_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, D_bank, Twater, Lambda, Tbank0(i), Tf, 1);

    [~,~,~, tssa_Tbank0_fit(i)] = SteadyStateSoln_fittedTss(Cf, U, S, D_bank, Twater, Lambda, Tbank0(i), Tf);
end

% get steady state conditions for Twater
load('YukonBaseCaseSS.mat');
load('Yukon_TLss_Twater.mat','eta_ss','Twater','eta_ts','M','E');
tss_Twater = eta_ts;
tssa_Twater = eta_ss./(M-E);
for i = 1:length(Twater)
    [Eent, ~, ~, eta_Twater_ss(i), Twater_Eents(i), Twater_Tws(i), Twater_Us(i), ...
        Twater_T0s(i), Twater_Cics(i), Twater_Chs(i), Twater_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U, S, D_bank, Twater(i), Lambda, Tbank0, Tf, 1);

    [~,~,~, tssa_Twater_fit(i)] = SteadyStateSoln_fittedTss(Cf, U, S, D_bank, Twater(i), Lambda, Tbank0, Tf);
end

% get steady state conditions for U
load('YukonBaseCaseSS.mat');
load('Yukon_TLss_U.mat','eta_ss','U','eta_ts','M','E');
tss_U = eta_ts;
tssa_U = eta_ss./(M-E);
for i = 1:length(U)
    [Eent, ~, ~, eta_U_ss(i), U_Eents(i), U_Tws(i), U_Us(i), ...
        U_T0s(i), U_Cics(i), U_Chs(i), U_Kts2s(i)] = ...
        SteadyStateSoln(Cf, U(i), S, D_bank, Twater, Lambda, Tbank0, Tf, 1);

    [~,~,~, tssa_U_fit(i)] = SteadyStateSoln_fittedTss(Cf, U(i), S, D_bank, Twater, Lambda, Tbank0, Tf);
end

%% Try to fit different regressions

Eents = [Dbank_Eents, Lambda_Eents, Tbank0_Eents, Twater_Eents, U_Eents]';
T0s = [Dbank_T0s, Lambda_T0s, Tbank0_T0s, Twater_T0s, U_T0s]';
Us = [Dbank_Us, Lambda_Us, Tbank0_Us, Twater_Us, U_Us]';
Chs = [Dbank_Chs, Lambda_Chs, Tbank0_Chs, Twater_Chs, U_Chs]';
Cics = [Dbank_Cics, Lambda_Cics, Tbank0_Cics, Twater_Cics, U_Cics]';
Tws = [Dbank_Tws, Lambda_Tws, Tbank0_Tws, Twater_Tws, U_Tws]';
Kts2s = [Dbank_Kts2s, Lambda_Kts2s, Tbank0_Kts2s, Twater_Kts2s, U_Kts2s]';
X1 = [Eents, abs(T0s), Us, Chs, Cics, Tws, Kts2s];

y = [tss_Dbank, tss_Lambda, tss_Tbank0, tss_Twater, tss_U]';
ymask = (~isinf(y) & ~isnan(y) & ~isinf(Eents) & ~isinf(T0s) & ~isinf(Us) ...
    & ~isinf(Chs) & ~isinf(Cics) & ~isinf(Tws) & ~isinf(Kts2s) & Eents>0 ...
    & T0s<0 & Tws>0);
y = y(ymask);
X1 = X1(ymask,:);
X2 = log10(abs(X1));
Tssfunc = fitlm(X2,log10(y));

save('FittedTssfunction.mat','Tssfunc');


%% comparison figure 11d

figure()
hold on
plot([0.01,100],[0.01,100],'k-');
plot(tss_Dbank/60/60, tssa_Dbank_fit/60/60,'o','LineWidth',2);
plot(tss_Lambda/60/60, tssa_Lambda_fit/60/60,'ks','LineWidth',2);
plot(tss_Tbank0/60/60, tssa_Tbank0_fit/60/60,'d','LineWidth',2);
plot(tss_Twater/60/60, tssa_Twater_fit/60/60,'>','LineWidth',2);
plot(tss_U/60/60, tssa_U_fit/60/60,'*','LineWidth',2);
xlabel('{\itt_s_s}');
ylabel('{\itt_s_s_a}');
set(gca,'FontSize',16,'YScale','log','XScale','log');
legend({'1:1 line','Vary {\itD_5_0}','Vary {\it\lambda_p}','Vary {\itT_{bank}}', ...
    'Vary {\itT_{water}}','Vary {\itU}'},'location','southeast');
box on
grid on
hold off
