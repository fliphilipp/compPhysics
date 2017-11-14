clear all;
close all;

% T / E_tot output from multiple MD simulations
solid = [489.98604, -806.20992;
490.97170, -806.15016;
491.96675, -806.05717;
492.98525, -805.94783;
493.94133, -805.88226;
494.95986, -805.81647;
495.95942, -805.72783;
496.96181, -805.63543;
497.92797, -805.58147;
498.95620, -805.49712;
499.94423, -805.40256;
500.92372, -805.32990;
501.93191, -805.25581;
502.92058, -805.16943;
503.95902, -805.07859;
504.92816, -805.00814;
505.91240, -804.95242;
506.87200, -804.86670;
507.94043, -804.72848;
508.92489, -804.64776;
509.89075, -804.61897];

fluid = [689.82340, -789.03003;
690.93559, -788.93400;
691.86218, -788.83357;
692.85667, -788.66970;
693.87984, -788.59014;
694.89975, -788.55241;
695.88103, -788.41978;
696.79565, -788.33073;
697.86977, -788.21197;
698.83720, -788.15122;
699.85247, -788.08793;
700.83904, -787.95794;
701.82991, -787.83477;
702.82874, -787.79467;
703.84887, -787.64322;
704.86329, -787.58728;
705.83564, -787.47751;
706.81277, -787.38276;
707.88622, -787.25960;
708.84338, -787.17674;
709.82887, -787.09849];

at20 = [14.99538, -840.77347;
15.99791, -840.70456;
16.99309, -840.63446;
17.99674, -840.57041;
18.99014, -840.50388;
19.99260, -840.43712;
20.99102, -840.36568;
21.99579, -840.29615;
22.99370, -840.22523;
23.99695, -840.15771;
24.99764, -840.09238;
25.99511, -840.02624;
26.98768, -839.95432;
27.99501, -839.88660;
28.99374, -839.82070;
29.99421, -839.75094;
30.99577, -839.68250;
32.00058, -839.62056;
32.99258, -839.54895;
33.99624, -839.47862;
34.99570, -839.41095];

xSolid = solid(:,1);
ySolid = solid(:,2);
pSolid = polyfit(xSolid,ySolid,1);
yfitSolid = polyval(pSolid,xSolid);

xFluid = fluid(:,1);
yFluid = fluid(:,2);
pFluid = polyfit(xFluid,yFluid,1);
yfitFluid = polyval(pFluid,xFluid);

fig1 = figure(1); set(fig1, 'Position', [100, 10, 1300, 700]);
subplot(1,2,1);
scatter(xSolid, ySolid);
hold on;
plot(xSolid, yfitSolid);
xlim([min(xSolid) max(xSolid)]);
ylim([min(ySolid) max(ySolid)]);
set(gca,'fontsize',18);
xlabel('Temperature $[^\circ C]$','interpreter','latex','fontsize',24);
ylabel('Total Energy $[eV]$','interpreter','latex','fontsize',24);
title('Total Energy vs Temperature around $500^\circ C$','interpreter','latex','fontsize',24)
fitdescriptionSolid = ['linear fit: ' '$\frac{\partial E_{tot}}{\partial T} \approx ' num2str(pSolid(1)) '\frac{eV}{K}$'];
lgd1 = legend('MD results', fitdescriptionSolid);
set(lgd1,'interpreter','latex','fontsize',20,'location','northwest')

subplot(1,2,2);
scatter(xFluid, yFluid);
hold on;
plot(xFluid, yfitFluid);
xlim([min(xFluid) max(xFluid)]);
ylim([min(yFluid) max(yFluid)]);
set(gca,'fontsize',18);
xlabel('Temperature $[^\circ C]$','interpreter','latex','fontsize',24);
ylabel('Total Energy $[eV]$','interpreter','latex','fontsize',24);
title('Total Energy vs Temperature around $700^\circ C$','interpreter','latex','fontsize',24)
fitdescriptionFluid = ['linear fit: ' '$\frac{\partial E_{tot}}{\partial T} \approx ' num2str(pFluid(1)) '\frac{eV}{K}$'];
lgd2 = legend('MD results', fitdescriptionFluid);
set(lgd2,'interpreter','latex','fontsize',20,'location','northwest')

setPDFsize;

saveas(fig1,'heatCapacityLinearFit.pdf')
saveas(fig1,'heatCapacityLinearFit.png')

x20 = at20(:,1);
y20 = at20(:,2);
p20 = polyfit(x20,y20,1);
yfit20 = polyval(p20,x20);

fig2 = figure(2); set(fig2, 'Position', [100, 10, 600, 500]);
scatter(x20, y20);
hold on;
plot(x20, yfit20);
xlim([min(x20) max(x20)]);
ylim([min(y20) max(y20)]);
set(gca,'fontsize',18);
xlabel('Temperature $[K]$','interpreter','latex','fontsize',24);
ylabel('Total Energy $[eV]$','interpreter','latex','fontsize',24);
title('Total Energy vs Temperature around $20^\circ C$','interpreter','latex','fontsize',24)
fitdescriptionFluid = ['linear fit: ' '$\frac{\partial E_{tot}}{\partial T} \approx ' num2str(p20(1)) '$'];
lgd3 = legend('MD results', fitdescriptionFluid);
set(lgd3,'interpreter','latex','fontsize',20,'location','northwest')

eliminateMargins;
setPDFsize;

saveas(fig1,'heatCapacityLinearFit20celsius.pdf')
saveas(fig1,'heatCapacityLinearFit20celsius.png')
