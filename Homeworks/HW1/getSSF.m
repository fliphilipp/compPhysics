%%
%clear all;
close all;

%%filenames
fnFluid = 'distances-fluid-700C.dat';
nBins = 2000;  % how many bins to use for the distances
splineP = 0.99999;  % smoothing parameter

%%
data = importdata(fnFluid);
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

%% get static structure factor according to Eq 58 in the lecture notes
qRange = 2.2:0.01:15;
ssf = zeros(length(qRange),1);
counter = 1;
for q = qRange
  forSSF = x.^2 .* (y - 1) .* sin(q * x) ./ (q * x); % tranform
  splineSSF = csaps(x,forSSF,splineP); % fit spline
  integral = diff(fnval(fnint(splineSSF),[0 cell_length])); % integrate up to cell length
  ssf(counter) = 1 + 4 * pi * 256 * integral;
  counter = counter + 1;
end

fig1 = figure(1);
plot(qRange, ssf, 'k-')
set(gca,'fontsize',14);
xlim([min(qRange) max(qRange)])
ylim([min(ssf) max(ssf)])
xlabel('$q$','fontsize',24,'interpreter','latex');
ylabel('$S(q)$ [arbitrary units]','interpreter','latex','fontsize',24);
title('Static structure factor, fluid Aluminum, $700^\circ C$','interpreter','latex','fontsize',24)
eliminateMargins;
setPDFsize
saveas(fig1,'static-structure-factor.pdf')