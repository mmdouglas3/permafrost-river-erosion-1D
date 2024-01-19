% plot results from steady-state models to produce figures 5 andf 6

%% bank grain size

clear
load('Yukon_TLss_Dbank.mat');
ind = find(Dbank>=D_bank,1);

figure()

subplot(3,4,1)
hold on
plot(Dbank,E,'k-','LineWidth',2);
plot(D_bank,E(ind),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itD_5_0} (m)');
ylabel('{\itE} (m/s)');
xlim([1e-5,0.02])
ylim([1e-10,1e-5])
xticks([1e-5,1e-4,1e-3,1e-2]);
yticks([10^(-10),1e-9,1e-8,1e-7,1e-6,1e-5])
box on
grid on
set(gca,'FontSize',14,'XScale','log','YScale','log');
text(-15,9,'a','FontSize',24,'Units','characters');
hold off

subplot(3,4,2)
hold on
plot(Dbank,q_star,'k-','LineWidth',2);
plot(D_bank,q_star(ind),'m*','LineWidth',2,'MarkerSize',20);
xlim([1e-5,0.02])
xticks([1e-5,1e-4,1e-3,1e-2]);
xlabel('{\itD_5_0} (m)');
ylabel('{\itq*}');
box on
grid on
set(gca,'FontSize',14,'XScale','log','YScale','log');
text(-15,9,'b','FontSize',24,'Units','characters');
hold off

subplot(3,4,3)
hold on
plot(Dbank,eta_ss./Dbank,'k-','LineWidth',2);
plot(D_bank,eta_ss(ind)./D_bank,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itD_5_0} (m)');
ylabel('{\it\eta*_s_s}');
xlim([1e-5,0.02])
xticks([1e-5,1e-4,1e-3,1e-2]);
box on
grid on
set(gca,'FontSize',14,'XScale','log','YScale','log');
text(-15,9,'c','FontSize',24,'Units','characters');
hold off

subplot(3,4,4)
M = sqrt(R*9.81*Dbank).*q_star;
tss = eta_ss./(M-E)/60/60;
hold on
plot(Dbank(~isnan(E)),eta_ts(~isnan(E))/60/60,'k-','LineWidth',2);
plot(D_bank,eta_ts(ind)/60/60,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itD_5_0} (m)');
ylabel('{\itt_s_s} (hr)');
xlim([1e-5,0.02])
xticks([1e-5,1e-4,1e-3,1e-2]);
box on
grid on
set(gca,'FontSize',14,'XScale','log','YScale','log');
text(-15,9,'d','FontSize',24,'Units','characters');
hold off

save('Dbank_SS.mat','Dbank','tss','eta_ts');

%% bank porosity

clear
load('Yukon_TLss_lambda.mat');
ind = find(Lp>=Lambda,1);

subplot(3,4,5)
hold on
plot(Lp,E,'k-','LineWidth',2);
plot(Lambda,E(ind),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\it\lambda_p}');
ylabel('{\itE} (m/s)');
ylim([1e-7,1e-5]);
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-15,9,'e','FontSize',24,'Units','characters');
hold off

subplot(3,4,6)
hold on
plot(Lp(1:end-1),q_star(1:end-1),'k-','LineWidth',2);
plot(Lambda,q_star(ind),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\it\lambda_p}');
ylabel('{\itq*}');
ylim([1e-4,1e-2]);
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-15,9,'f','FontSize',24,'Units','characters');
hold off

subplot(3,4,7)
hold on
plot(Lp,eta_ss./D_bank,'k-','LineWidth',2);
plot(Lambda,eta_ss(ind)/D_bank,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\it\lambda_p}');
ylabel('{\it\eta*_s_s}');
ylim([1,1e3]);
yticks([1,1e1,1e2,1e3]);
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-15,9,'g','FontSize',24,'Units','characters');
hold off

subplot(3,4,8)
M = sqrt(R*9.81*D_bank).*q_star;
tss = eta_ss./(M-E)/60/60;
hold on
plot(Lp(2:end-1),eta_ts(2:end-1)/60/60,'k-','LineWidth',2);
plot(Lambda,eta_ts(ind)/60/60,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\it\lambda_p}');
ylabel('{\itt_s_s} (hr)');
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-15,9,'h','FontSize',24,'Units','characters');
hold off

save('Lambda_SS.mat','Lp','tss','eta_ts');

%% bank temperature

clear
load('Yukon_TLss_Tbank.mat');

subplot(3,4,9)
M = sqrt(R*9.81*D_bank).*q_star;
hold on
plot(Tbank0,E,'k-','LineWidth',2);
plot(-1,E(Tbank0==-1),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_0} (\circC)');
ylabel('{\itE} (m/s)');
ylim([0,1e-6]);
box on
grid on
set(gca,'FontSize',14);
text(-14,9,'i','FontSize',24,'Units','characters');
hold off

subplot(3,4,10)
hold on
plot(Tbank0,q_star,'k-','LineWidth',2);
plot(-1,q_star(Tbank0==-1),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_0} (\circC)');
ylabel('{\itq*}');
ylim([0,1e-2])
box on
grid on
set(gca,'FontSize',14);
text(-14,9,'j','FontSize',24,'Units','characters');
hold off

subplot(3,4,11)
hold on
plot(Tbank0,eta_ss./D_bank,'k-','LineWidth',2);
plot(-1,eta_ss(Tbank0==-1)/D_bank,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_0} (\circC)');
ylabel('{\it\eta*_s_s}');
ylim([0,100])
box on
grid on
set(gca,'FontSize',14);
text(-14,9,'k','FontSize',24,'Units','characters');
hold off

subplot(3,4,12)
M = sqrt(R*9.81*D_bank).*q_star;
tss = eta_ss./(M-E)/60/60;
hold on
plot(Tbank0(~isnan(E)),eta_ts(~isnan(E))/60/60,'k-','LineWidth',2);
plot(-1,eta_ts(Tbank0==-1)/60/60,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_0} (\circC)');
ylabel('{\itt_s_s} (hr)');
box on
grid on
set(gca,'FontSize',14);
text(-14,9,'l','FontSize',24,'Units','characters');
hold off

save('Tbank0_SS.mat','Tbank0','tss','eta_ts');

%% water temperature

clear
load('Yukon_TLss_Twater.mat');
ind = find(Twater>=10,1);

figure()

subplot(2,4,1)
hold on
plot(Twater,E,'k-','LineWidth',2);
plot(10,E(ind),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_{w}} (\circC)');
ylabel('{\itE} (m/s)');
ylim([1e-7,1e-6])
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-14,11,'a','FontSize',24,'Units','characters');
hold off

subplot(2,4,2)
hold on
plot(Twater,q_star,'k-','LineWidth',2);
plot(10,q_star(ind),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_{w}} (\circC)');
ylabel('{\itq*}');
ylim([1e-4,1e-2]);
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-14,11,'b','FontSize',24,'Units','characters');
hold off

subplot(2,4,3)
hold on
plot(Twater,eta_ss./D_bank,'k-','LineWidth',2);
plot(10,eta_ss(ind)/D_bank,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_{w}} (\circC)');
ylabel('{\it\eta*_s_s}');
ylim([0,120]);
box on
grid on
set(gca,'FontSize',14);
text(-14,11,'c','FontSize',24,'Units','characters');
hold off

subplot(2,4,4)
M = sqrt(R*9.81*D_bank).*q_star;
tss = eta_ss./(M-E)/60/60;
hold on
plot(Twater(~isnan(E)),eta_ts(~isnan(E))/60/60,'k-','LineWidth',2);
plot(10,eta_ts(ind)/60/60,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itT_{w}} (\circC)');
ylabel('{\itt_s_s} (hr)');
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-14,11,'d','FontSize',24,'Units','characters');
hold off

save('Twater_SS.mat','Twater','tss','eta_ts');

%% water velocity

clear
load('Yukon_TLss_U.mat');

subplot(2,4,5)
hold on
plot(U,E,'k-','LineWidth',2);
plot(1,E(U==1),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itU} (m/s)');
ylabel('{\itE} (m/s)');
ylim([1e-7,1e-5])
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-14,11,'e','FontSize',24,'Units','characters');
hold off

subplot(2,4,6)
hold on
plot(U,q_star,'k-','LineWidth',2);
plot(1,q_star(U==1),'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itU} (m/s)');
ylabel('{\itq*}');
ylim([1e-3,3e-3])
box on
grid on
set(gca,'FontSize',14);
text(-14,11,'f','FontSize',24,'Units','characters');
hold off

subplot(2,4,7)
hold on
plot(U,eta_ss./D_bank,'k-','LineWidth',2);
plot(1,eta_ss(U==1)/D_bank,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itU} (m/s)');
ylabel('{\it\eta*_s_s}');
box on
grid on
set(gca,'FontSize',14);
text(-14,11,'g','FontSize',24,'Units','characters');
hold off

subplot(2,4,8)
M = sqrt(R*9.81*D_bank).*q_star;
tss = eta_ss./(M-E)/60/60;
hold on
plot(U(~isnan(E)),eta_ts(~isnan(E))/60/60,'k-','LineWidth',2);
plot(1,eta_ts(U==1)/60/60,'m*','LineWidth',2,'MarkerSize',20);
xlabel('{\itU} (m/s)');
ylabel('{\itt_s_s} (hr)');
box on
grid on
set(gca,'FontSize',14,'YScale','log');
text(-14,11,'h','FontSize',24,'Units','characters');
hold off

save('U_SS.mat','U','tss','eta_ts');

