% make Figure 8 with results for Yukon River summer flood

clear

int = 1e5;
load('YukonHydrographSF.mat','Qw','Tw');
load('YukonHydrograph_modernSF.mat','EP','MP','icP','tsP','doy','dt');
Ebase = EP(1:int:end);
Ebase(Ebase<0) = 0;
MP(find(MP>0,1):find(MP>0,1,'last')) = movmean(MP(find(MP>0,1):find(MP>0,1,'last')),1e6);
Mbase = MP(1:int:end);
Mbase(Mbase<0) = 0;
ebase = min([Ebase,Mbase],[],2);
eta_base = tsP(1:int:end) - icP(1:int:end);
doy = transpose(doy(2:int:end));
Qw = Qw(2:int:end)';
Tw = Tw(2:int:end)';
Tw(Tw<0) = 0;

figure('Renderer', 'painters', 'Position', [0 0 1100 1000])

subplot(5,1,1)
hold on
plot(doy,Qw,'k','LineWidth',1.5);
ylabel('{\itQ_w} (m^3/s)');
ylim([0,11e3]);
xlim([110,310]);
set(gca,'FontSize',16);
text(1,3,'a','Units','characters','FontSize',24)
box on
hold off

subplot(5,1,2)
hold on
plot(doy,Tw,'k','LineWidth',1.5);
ylabel('{\itT_w} (\circC)');
xlim([110,310]);
set(gca,'FontSize',16);
text(1,3,'b','Units','characters','FontSize',24)
box on
hold off

subplot(5,1,3)
hold on
plot(doy,Mbase,'m','LineWidth',1);
plot(doy,Ebase,'b','LineWidth',1);
xlim([110,310]);
ylabel('{\itE} (m/s)');
legend({'{\itE_{thaw}}','{\itE_{ent}}','{\itE}'});
set(gca,'FontSize',16);
text(1,3,'c','Units','characters','FontSize',24)
box on
hold off

subplot(5,1,4)
hold on
plot(doy,tsP(1) - tsP(1:int:end),'k-','LineWidth',2);
ylabel('Dist. eroded (m)');
xlim([110,310]);
set(gca,'FontSize',16);
text(1,3,'d','Units','characters','FontSize',24)
box on
hold off

subplot(5,1,5)
hold on
plot(doy,eta_base,'k-','LineWidth',2);
ylabel('\eta (m)');
xlim([110,310]);
set(gca,'FontSize',16);
text(1,3,'e','Units','characters','FontSize',24)
box on
hold off
