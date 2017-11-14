%%
%clear all;
close all;

%%filenames
fnSolid = 'distances-solid-20C.dat';
fnFluid = 'distances-fluid-700C.dat';
nBins = 2000;  % how many bins to use for the distances
splineP = 0.99999;  % smoothing parameter

%% solid
data = importdata(fnSolid);
disp('Solid')
nMeasurements = data(end);  % wrote this to last position in C
cell_length = data(end-1);  % wrote this to second last position in C
disp(['solid measurements: ' num2str(nMeasurements)])
disp(['solid cell length: ' num2str(cell_length)])
dist = data(1:end-2);  % the actual data of all distances
disp(['solid number of values: ' num2str(length(dist))])
density = 256 * nMeasurements / cell_length^3;
disp(['solid density: ' num2str(density)])

%% bin into histogram counts
[N,r] = histcounts(dist,nBins);
dr = r(2) - r(1);
nmissing = min(r) / dr;
addr = (min(r)-nmissing*dr):dr:(min(r)-dr);
r = [addr, r];
N = [zeros(1,length(addr)) N];
scaling = 4 / 3 * pi * ((r + dr).^3 - r.^3) * density * 256;
scaling(end) = [];
middle = r + dr/2;
x = middle(1:end-1);
y = N ./ scaling;
smoothspline = csaps(x,y,splineP);
[minval,minsite] = fnmin(smoothspline,[3 4]);
[peak1val,peak1site] = fnmin(fncmb(smoothspline,-1),[0 3]);
[peak2val,peak2site] = fnmin(fncmb(smoothspline,-1),[3.8 4.3]);
peak1val = - peak1val;
peak2val = - peak2val;
disp(['solid first peak at ' num2str(peak1site)])
disp(['solid second peak at ' num2str(peak2site)])

%% calculate coordination number
forspline = y .* x.^2; % multiply by r^2 for coordination number integration
spliner2 = csaps(x,forspline,splineP);
coord = density / nMeasurements * 4 * pi * diff(fnval(fnint(spliner2),[0 minsite]));
disp(['Coordination number: ' num2str(coord,8)])

%% plot
fig1 = figure(1); set(fig1, 'Position', [100, 10, 1000, 700]);
subplot(2,1,1)
scatter(x, y, [], 'k')
hold on
plot([0 30], [1 1], 'b--', 'linewidth', 3)
fnplt(smoothspline,'y-')
xlim([0 16])
ylim([0 9])
scatter(minsite,minval,200,'red','filled')
scatter(peak1site,peak1val,200,'green','filled')
scatter(peak2site,peak2val,200,'green','filled')
descrip = {['$r_m = ' num2str(minsite) '$'],['$I(r_m) = ' num2str(coord) '$']};
text(4,7,descrip,'interpreter','latex','fontsize',28)
text(peak1site,peak1val+1,[num2str(peak1site,3) 'Å'],'fontsize',20,'HorizontalAlignment','center')
text(peak2site,peak2val+1,[num2str(peak2site,3) 'Å'],'fontsize',20,'HorizontalAlignment','center')
set(gca,'fontsize',18);
xlabel('distance [Å]','fontsize',20);
ylabel('$g(r)$','interpreter','latex','fontsize',30);
title('solid Aluminum, $20^\circ C$','interpreter','latex','fontsize',30)
lgd1 = legend('MD results', '$g(r) = 1$', 'spline fit','first minimum of $g(r)$','face-centered cubic first peaks');
set(lgd1,'interpreter','latex','fontsize',20,'location','northeast')

%% fluid
data = importdata(fnFluid);
disp('')
disp('Fluid')
nMeasurements = data(end);  % wrote this to last position in C
cell_length = data(end-1);  % wrote this to second last position in C
disp(['fluid measurements: ' num2str(nMeasurements)])
disp(['fluid cell length: ' num2str(cell_length)])
dist = data(1:end-2);  % the actual data of all distances
dist(dist > cell_length) = [];
disp(['fluid number of values: ' num2str(length(dist))])
density = 256 * nMeasurements / cell_length^3;
disp(['fluid density: ' num2str(density)])

%% bin into histogram counts
[N,r] = histcounts(dist,nBins);
dr = r(2) - r(1);
nmissing = min(r) / dr;
addr = (min(r)-nmissing*dr):dr:(min(r)-dr);
r = [addr, r];
N = [zeros(1,length(addr)) N];
scaling = 4 / 3 * pi * ((r + dr).^3 - r.^3) * density * 256;
scaling(end) = [];
middle = r + dr/2;
x = middle(1:end-1);
y = N ./ scaling;
smoothspline = csaps(x,y,0.99999);
[minval,minsite] = fnmin(smoothspline,[3 5]);

%% calculate coordination number
forspline = y .* x.^2; % multiply by r^2 for coordination number integration
spliner2 = csaps(x,forspline,splineP);
coord = density / nMeasurements * 4 * pi * diff(fnval(fnint(spliner2),[0 minsite]));
disp(['Coordination number: ' num2str(coord,8)])

%% plot
subplot(2,1,2)
scatter(x, y, [], 'k')
hold on
xlim([0,16])
ylim([0 4])
plot([0 30], [1 1], 'b--', 'linewidth', 3)
fnplt(smoothspline,'y-')
scatter(minsite,minval,200,'red','filled')
descrip = {['$r_m = ' num2str(minsite) '$'],['$I(r_m) = ' num2str(coord) '$']};
text(4,3,descrip,'interpreter','latex','fontsize',28)
set(gca,'fontsize',18);
xlabel('distance [Å]','fontsize',20);
ylabel('$g(r)$','interpreter','latex','fontsize',30);
title('fluid Aluminum, $700^\circ C$','interpreter','latex','fontsize',30)
lgd1 = legend('MD results', '$g(r) = 1$', 'spline fit','first minimum of $g(r)$');
set(lgd1,'interpreter','latex','fontsize',20,'location','northeast')

saveas(fig1,'radial-dist-func.png')

