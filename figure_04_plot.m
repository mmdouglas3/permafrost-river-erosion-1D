% plot figure 4

clear

tmult = 1e4;
xmult = 1;

load('Yukon_sand_lambda50.mat','dx','bankdepth','dt','totalt','icP','tsP','EP');
ic_s50 = icP;
ts_s50 = tsP;
L_s50 = ts_s50 - ic_s50;
Tss_s50 = L_s50(end)/EP(end)/60/60;
tt = transpose(dt:dt*tmult:totalt)/60/60;
xx = transpose(0:dx*xmult:bankdepth+dx);
clear('dx','bankdepth','dt','totalt','icP','tsP','EP');

load('Yukon_sand_lambda22.mat','icP','tsP','EP');
ic_s22 = icP;
ts_s22 = tsP;
L_s22 = ts_s22 - ic_s22;
Tss_s22 = L_s22(end)/EP(end)/60/60;

load('Yukon_sand_lambda80.mat','icP','tsP','EP');
ic_s80 = icP;
ts_s80 = tsP;
L_s80 = ts_s80 - ic_s80;
Tss_s80 = L_s80(end)/EP(end)/60/60;

load('Yukon_silt_lambda22.mat','icP','tsP','EP');
ic_l22 = icP;
ts_l22 = tsP;
L_l22 = ts_l22 - ic_l22;
Tss_l22 = L_l22(end)/EP(end)/60/60;

load('Yukon_gravel_lambda22.mat','icP','tsP','EP');
ic_g22 = icP;
ts_g22 = tsP;
L_g22 = ts_g22 - ic_g22;
Tss_g22 = L_g22(end)/EP(end)/60/60;
clear('dx','bankdepth','dt','totalt','icP','tsP','EP');

load('Yukon_gravel_lambda50.mat','icP','tsP','EP');
ic_g50 = icP;
ts_g50 = tsP;
L_g50 = ts_g50 - ic_g50;
Tss_g50 = L_g50(end)/EP(end)/60/60;
clear('dx','bankdepth','dt','totalt','icP','tsP','EP');

load('Yukon_gravel_lambda80.mat','icP','tsP','EP');
ic_g80 = icP;
ts_g80 = tsP;
L_g80 = ts_g80 - ic_g80;
Tss_g80 = L_g80(end)/EP(end)/60/60;
clear('dx','bankdepth','dt','totalt','icP','tsP','EP');

figure('Renderer', 'painters', 'Position', [10 10 700 1100])
legtext = {'Permafrost','Thawed sediment'};
c1 = [0 0.5 1];
c2 = [0.7 0.7 0.7];
xlimits = [0,48];
ylimits = [4.6,5];
yt = 0:0.1:0.4;

subplot(7,1,1);
hold on
h = area(tt, [ic_s22(1:tmult:end), L_s22(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'a','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 1 mm, {\it\lambda_p} = 0.22','Units','characters','FontSize',11)
box on
hold off

subplot(7,1,2);
hold on
h = area(tt, [ic_s50(1:tmult:end), L_s50(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
legend(legtext,'location','southeast','orientation','horizontal');
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'b','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 1 mm, {\it\lambda_p} = 0.5','Units','characters','FontSize',11)
box on
hold off

subplot(7,1,3);
hold on
h = area(tt, [ic_s80(1:tmult:end), L_s80(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'c','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 1 mm, {\it\lambda_p} = 0.8','Units','characters','FontSize',11)
box on
hold off

subplot(7,1,4);
hold on
h = area(tt, [ic_l22(1:tmult:end), L_l22(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'d','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 50 \mum, {\it\lambda_p} = 0.22','Units','characters','FontSize',11)
box on
hold off

subplot(7,1,5);
hold on
h = area(tt, [ic_g22(1:tmult:end), L_g22(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'e','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 1 cm, {\it\lambda_p} = 0.22','Units','characters','FontSize',11)
box on
hold off

subplot(7,1,6);
hold on
h = area(tt, [ic_g50(1:tmult:end), L_g50(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'f','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 1 cm, {\it\lambda_p} = 0.5','Units','characters','FontSize',11)
box on
hold off

subplot(7,1,7);
hold on
h = area(tt, [ic_g80(1:tmult:end), L_g80(1:tmult:end)]);
h(1).FaceColor = c1;
h(2).FaceColor = c2;
ylabel('{\itX} (m)');
xlabel('{\itt} (hr)')
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',10);
text(-10,3,'g','Units','characters','FontSize',16)
text(1,1,'{\itD_{50}} = 1 cm, {\it\lambda_p} = 0.8','Units','characters','FontSize',11)
box on
hold off


