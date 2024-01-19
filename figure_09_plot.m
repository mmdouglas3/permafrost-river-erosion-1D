% plot figure 9 with Yukon hydrograph bank failure results

clear

int = 1e3;
load('YukonHydrographBF.mat','Qw','Tw');
load('YukonHydrograph_BankFail.mat','EP','MP','icP','tsP','doy','dt');
Ebase = EP(1:int:end);
Ebase(Ebase<0) = 0;
MP(find(MP>0,1):find(MP>0,1,'last')) = movmean(MP(find(MP>0,1):find(MP>0,1,'last')),1e5);
Mbase = MP(1:int:end);
Mbase(Mbase<0) = 0;
eroded = 50-tsP(1:int:end);
eta_base = tsP - icP;
eta_daily = movmean(eta_base,24*60*6*5);
doy2 = transpose(doy(2:end));
doy = transpose(doy(2:int:end));
Qw = Qw(2:int:end)';
Tw = Tw(2:int:end)';
Tw(Tw<0) = 0;

figure('Renderer', 'painters', 'Position', [0 0 1100 600])

subplot(3,1,1)
hold on
plot(doy,Mbase,'m','LineWidth',1);
plot(doy,Ebase,'b','LineWidth',1);
xlim([110,310]);
ylabel('{\itE} (m/s)');
legend({'{\itE_{thaw}}','{\itE_{ent}}','{\itE}'});
set(gca,'FontSize',16);
text(1,8,'a','Units','characters','FontSize',24)
box on
hold off

subplot(3,1,2)
hold on
plot(doy,eroded,'k-','LineWidth',2);
ylabel('Dist. eroded (m)');
xlim([110,310]);
set(gca,'FontSize',16);
text(1,8,'b','Units','characters','FontSize',24)
box on
hold off

subplot(3,1,3)
hold on
plot(doy2,eta_base,'k-','LineWidth',2);
plot(doy2,eta_daily,'c-','LineWidth',2);
legend({'instantaneous','5-day mean'},'location','south','Orientation','horizontal')
ylabel('\eta (m)');
xlim([110,310]);
set(gca,'FontSize',16);
text(1,8,'c','Units','characters','FontSize',24)
box on
hold off
