% plot bank failure results

clear

bankshift = 0;
timeint = 1e4;

%% load 1 cm D50

load('Yukon_nofail_D1cm.mat','icP','tsP');
ic_nfcm = icP(1:timeint:end)-bankshift;
ts_nfcm = tsP(1:timeint:end)-bankshift;

load('Yukon_1cmFail.mat','icP','tsP');
ic_1cm = icP(1:timeint:end)-bankshift;
ts_1cm = tsP(1:timeint:end)-bankshift;

load('Yukon_3cmFail.mat','icP','tsP');
ic_3cm = icP(1:timeint:end)-bankshift;
ts_3cm = tsP(1:timeint:end)-bankshift;

load('Yukon_5cmFail.mat','icP','tsP');
ic_5cm = icP(1:timeint:end)-bankshift;
ts_5cm = tsP(1:timeint:end)-bankshift;

load('Yukon_10cmFail.mat','icP','tsP');
ic_10cm = icP(1:timeint:end)-bankshift;
ts_10cm = tsP(1:timeint:end)-bankshift;

load('Yukon_15cmFail.mat','icP','tsP');
ic_15cm = icP(1:timeint:end)-bankshift;
ts_15cm = tsP(1:timeint:end)-bankshift;

load('Yukon_25cmFail.mat','icP','tsP','dt','totalt');
ic_25cm = icP(1:timeint:end)-bankshift;
ts_25cm = tsP(1:timeint:end)-bankshift;
tt = dt:(dt*timeint):totalt;
tt = tt/24/60/60;
clear icP tsP dt totalt

yFail_cm = [1;3;5;10;15;25;100*(ts_nfcm(end) - ic_nfcm(end))];
eta_ratio_cm = yFail_cm/(ts_nfcm(end) - ic_nfcm(end))/100;
Erate_cm = (ts_nfcm(1) - [ts_1cm(end); ts_3cm(end); ts_5cm(end); ts_10cm(end); ...
    ts_15cm(end); ts_25cm(end); ts_nfcm(end)]) / max(tt(:));

%% load 1 mm D50

load('Yukon_nofail_D1mm.mat','icP','tsP');
ic_nfmm = icP(1:timeint:end)-bankshift;
ts_nfmm = tsP(1:timeint:end)-bankshift;

load('Yukon_1mmFail.mat','icP','tsP');
ic_1mm = icP(1:timeint:end)-bankshift;
ts_1mm = tsP(1:timeint:end)-bankshift;

load('Yukon_30mmFail.mat','icP','tsP');
ic_30mm = icP(1:timeint:end)-bankshift;
ts_30mm = tsP(1:timeint:end)-bankshift;

load('Yukon_50mmFail.mat','icP','tsP');
ic_50mm = icP(1:timeint:end)-bankshift;
ts_50mm = tsP(1:timeint:end)-bankshift;

load('Yukon_100mmFail.mat','icP','tsP');
ic_100mm = icP(1:timeint:end)-bankshift;
ts_100mm = tsP(1:timeint:end)-bankshift;

load('Yukon_150mmFail.mat','icP','tsP');
ic_150mm = icP(1:timeint:end)-bankshift;
ts_150mm = tsP(1:timeint:end)-bankshift;

load('Yukon_250mmFail.mat','icP','tsP');
ic_250mm = icP(1:timeint:end)-bankshift;
ts_250mm = tsP(1:timeint:end)-bankshift;
clear icP tsP

yFail_mm = [10;30;50;100;150;250;1000*(ts_nfmm(end) - ic_nfmm(end))];
eta_ratio_mm = yFail_mm/(ts_nfmm(end) - ic_nfmm(end))/1000;
Erate_mm = (ts_nfmm(1) - [ts_1mm(end); ts_30mm(end); ts_50mm(end); ts_100mm(end); ...
    ts_150mm(end); ts_250mm(end); ts_nfmm(end)]) / max(tt(:));

%% plot results

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(2,2,1)
hold on
plot(tt,ts_nfcm,'k-','LineWidth',3)
plot(tt,ts_1cm,'k-','LineWidth',2)
plot(tt,ts_3cm,'k-.','LineWidth',2)
plot(tt,ts_5cm,'k:','LineWidth',2)
plot(tt,ts_10cm,'k--','LineWidth',2)
legend({'no failures','1 cm','3 cm','5 cm','10 cm'},'location','southwest')
xlabel('Time (days)')
ylabel('{\itX} (m)')
set(gca,'FontSize',16)
text(0.6,1.5,'{\itD_5_0} = 1 cm','FontSize',16)
text(-0.2,15,'a','FontSize',24)
ylim([0,15]);
box on
hold off

subplot(2,2,3)
hold on
plot(tt,ts_nfmm,'k-','LineWidth',3)
plot(tt,ts_1mm,'k-','LineWidth',2)
plot(tt,ts_30mm,'k-.','LineWidth',2)
plot(tt,ts_50mm,'k:','LineWidth',2)
plot(tt,ts_100mm,'k--','LineWidth',2)
legend({'no failures','1 mm','3 cm','5 cm','10 cm'},'location','southwest')
xlabel('Time (days)')
ylabel('{\itX} (m)')
set(gca,'FontSize',16)
text(0.6,1.5,'{\itD_5_0} = 1 mm','FontSize',16)
text(-0.2,15,'c','FontSize',24)
ylim([0,15]);
box on
hold off

subplot(2,2,4)
hold on
plot([1,1],[1e-2,1e2],'k-','LineWidth',2);
plot(eta_ratio_mm,Erate_mm,'ks','MarkerFaceColor','w','MarkerSize',10,'LineWidth',2)
xlabel('{\it\eta_{fail} / \eta_{ss}}')
ylabel('{\itE} (m/day)')
set(gca,'FontSize',16,'YScale','log')
text(1.75,50,'{\itD_5_0} = 1 mm','FontSize',16)
text(-0.5,1e2,'d','FontSize',24)
ylim([1e-2,1e2]);
box on
hold off


%% load data for 1C Tw

load('Yukon_nofail_Tw1C.mat','icP','tsP');
eta_nf_Tw1 = tsP(end) - icP(end);
yFail_cm = [0.1;0.2;0.3;0.4;0.5;1;100*eta_nf_Tw1];
eta_ratio_Tw1 = yFail_cm/eta_nf_Tw1/100;
Erate_Tw1 = zeros(size(yFail_cm));
Erate_Tw1(7) = (tsP(1) - tsP(end))/ max(tt(:));
clear icP tsP

load('Yukon_1cmFail_Tw1C.mat','tsP');
Erate_Tw1(1) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_3cmFail_Tw1C.mat','tsP');
Erate_Tw1(2) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_5cmFail_Tw1C.mat','tsP');
Erate_Tw1(3) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_10cmFail_Tw1C.mat','tsP');
Erate_Tw1(4) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_15cmFail_Tw1C.mat','tsP');
Erate_Tw1(5) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_25cmFail_Tw1C.mat','tsP');
Erate_Tw1(6) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

%% load data for 5C Tw

load('Yukon_nofail_Tw5C.mat','icP','tsP');
eta_nf_Tw5 = tsP(end) - icP(end);
yFail_cm = [0.1;0.2;0.3;0.4;0.5;1;100*eta_nf_Tw5];
eta_ratio_Tw5 = yFail_cm/eta_nf_Tw5/100;
Erate_Tw5 = zeros(size(yFail_cm));
Erate_Tw5(7) = (tsP(1) - tsP(end))/ max(tt(:));
clear icP tsP

load('Yukon_1cmFail_Tw5C.mat','tsP');
Erate_Tw5(1) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_3cmFail_Tw5C.mat','tsP');
Erate_Tw5(2) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_5cmFail_Tw5C.mat','tsP');
Erate_Tw5(3) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_10cmFail_Tw5C.mat','tsP');
Erate_Tw5(4) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_15cmFail_Tw5C.mat','tsP');
Erate_Tw5(5) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_25cmFail_Tw5C.mat','tsP');
Erate_Tw5(6) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP


%% load data for 15C Tw

load('Yukon_nofail_Tw15C.mat','icP','tsP');
eta_nf_Tw15 = tsP(end) - icP(end);
yFail_cm = [0.1;0.2;0.3;0.4;0.5;1;100*eta_nf_Tw15];
eta_ratio_Tw15 = yFail_cm/eta_nf_Tw15/100;
Erate_Tw15 = zeros(size(yFail_cm));
Erate_Tw15(7) = (tsP(1) - tsP(end))/ max(tt(:));
clear icP tsP

load('Yukon_1cmFail_Tw15C.mat','tsP');
Erate_Tw15(1) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_3cmFail_Tw15C.mat','tsP');
Erate_Tw15(2) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_5cmFail_Tw15C.mat','tsP');
Erate_Tw15(3) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_10cmFail_Tw15C.mat','tsP');
Erate_Tw15(4) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_15cmFail_Tw15C.mat','tsP');
Erate_Tw15(5) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

load('Yukon_25cmFail_Tw15C.mat','tsP');
Erate_Tw15(6) = (tsP(1) - tsP(end))/ max(tt(:));
clear tsP

%% plot temperature dependence

M = [Erate_Tw1, Erate_Tw5, Erate_cm, Erate_Tw15];
style = {'ks-','ko--','kd:','k>-.','','','kx-'};
subplot(2,2,2)
hold on
for i = [1:4,7]
    plot([1,5,10,15],M(i,:),style{i},'LineWidth',2,'MarkerFaceColor','w','MarkerSize',5)
end
legend({'\eta_{fail} = 1 cm','\eta_{fail} = 3 cm','\eta_{fail} = 5 cm', ...
    '\eta_{fail} = 10 cm','no failures'});
xlabel('{\itT_w} (\circC)')
ylabel('{\itE} (m/day)')
set(gca,'FontSize',16)
text(10,9,'{\itD_5_0} = 1 cm','FontSize',16)
text(-3,10,'b','FontSize',24)
ylim([0,10]);
box on
hold off

