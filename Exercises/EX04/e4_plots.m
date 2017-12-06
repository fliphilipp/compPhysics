%% data import
clear all;
close all;

data_time = dlmread('time.data');

data_Ax = dlmread('trajAx.data');
data_Av = dlmread('trajAv.data') .* 1000;

data_Bx = dlmread('trajBx.data');
data_Bv = dlmread('trajBv.data') .* 1000;

data_A = importdata('trajA_Average.data');
data_A(:,3:4) = data_A(:,3:4) .* 1000;
data_B = importdata('trajB_Average.data');
data_B(:,3:4) = data_B(:,3:4) .* 1000;

ntraj = size(data_Ax, 2);
nmean = 10000;

xA = importdata('xdistA.data');
vA = importdata('vdistA.data') .* 1000;
xB = importdata('xdistB.data');
vB = importdata('vdistB.data') .* 1000;
nA = size(xA,1);
nB = size(xB,1);

pspec = importdata('pspec.data');
f_0 = 0.003; % frequency in MHz

%% task 2 plot
fig1 = figure(1); set(fig1, 'Position', [100, 10, 1200, 800],'DefaultAxesPosition', [0.03, 0.03, 0.97, 0.97]);

% positions A
subplot(2,2,1)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Ax(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_A(:,1) + data_A(:,2), 'b:','linewidth',2);
plot(data_time, data_A(:,1) - data_A(:,2), 'b:','linewidth',2)
h2 = plot(data_time, data_A(:,1),'r-','linewidth', 3);

disp(' ')
disp('get max/min values:')

lasti = 1;
for ii = 1:2
  [m, i] = max(data_A(lasti:end,1));
  im = data_time(lasti + i);
  disp(['max XA at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;

  [m, i] = min(data_A(lasti:end,1));
  im = data_time(lasti + i);
  disp(['min XA at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;
end
h4 = scatter(data_time(end),data_A(end,1),70,'g','filled');

set(gca, 'fontsize', 14)
title('Positions for relaxation time \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend([h1 h2 h3 h4], {'sample trajectories', ...
  ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation',...
  'distribution measurements'});
set(lgd,'fontsize',14)

% velocities A
subplot(2,2,2)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Av(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_A(:,3) + data_A(:,4), 'b:','linewidth',2);
plot(data_time, data_A(:,3) - data_A(:,4), 'b:','linewidth',2)
h2 = plot(data_time, data_A(:,3),'r-','linewidth', 3);

disp(' ')
lasti = 1;
for ii = 1:2
  [m, i] = max(data_A(lasti:end,3));
  im = data_time(lasti + i);
  disp(['max XV at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;

  [m, i] = min(data_A(lasti:end,3));
  im = data_time(lasti + i);
  disp(['min XV at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;
end
h4 = scatter(data_time(end),data_A(end,3),70,'g','filled');

set(gca, 'fontsize', 14)
title('Velocities for relaxation time \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend([h1 h2 h3], {'sample trajectories', ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation'});
set(lgd,'fontsize',14)
lgd = legend([h1 h2 h3 h4], {'sample trajectories',...
  ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation',...
  'distribution measurements'});

% positions B
subplot(2,2,3)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Bx(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_B(:,1) + data_B(:,2), 'b:','linewidth',2);
plot(data_time, data_B(:,1) - data_B(:,2), 'b:','linewidth',2)
h2 = plot(data_time, data_B(:,1),'r-','linewidth', 3);

disp(' ')
lasti = 1;
for ii = 1:4
  [m, i] = max(data_B(lasti:end,1));
  im = data_time(lasti + i);
  disp(['max XB at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;

  [m, i] = min(data_B(lasti:end,1));
  im = data_time(lasti + i);
  disp(['min XB at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;
end
h4 = scatter(data_time(end),data_B(end,1),70,'g','filled');

set(gca, 'fontsize', 14)
title('Positions for relaxation time \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend([h1 h2 h3 h4], {'sample trajectories',...
  ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation',...
  'distribution measurements'});
set(lgd,'fontsize',14)

% velocities B
subplot(2,2,4)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Bv(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_B(:,3) + data_B(:,4), 'b:','linewidth',2);
plot(data_time, data_B(:,3) - data_B(:,4), 'b:','linewidth',2)
h2 = plot(data_time, data_B(:,3),'r-','linewidth', 3);

disp(' ')
lasti = 1;
for ii = 1:4
  [m, i] = max(data_B(lasti:end,3));
  im = data_time(lasti + i);
  disp(['max XV at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;

  [m, i] = min(data_B(lasti:end,3));
  im = data_time(lasti + i);
  disp(['min XV at ' num2str(im)])
  scatter(im,m,70,'g','filled')
  lasti = lasti + i;
end
h4 = scatter(data_time(end),data_B(end,3),70,'g','filled');

set(gca, 'fontsize', 14)
title('Velocities for relaxation time \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend([h1 h2 h3 h4], {'sample trajectories',...
  ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation',...
  'distribution measurements'});
set(lgd,'fontsize',14)

saveas(fig1, 'task2.png')
setPDFsize;
saveas(fig1, 'task2.pdf')

%% Task 3 plotting
fig2 = figure(2); set(fig2, 'Position', [0, 0, 1400, 800],'DefaultAxesPosition', [0.03, 0.03, 0.97, 0.97]);
cols = get(gca,'colororder');
cols = [cols; cols];

txa = string({'$t=32.4\mu s$', '$t=231.5\mu s$', '$t=434.1\mu s$', '$t=645.3\mu s$', '$t=2000\mu s$'});
tva = string({'$t=0.1\mu s$', '$t=94.7\mu s$', '$t=293.6\mu s$', '$t=499.6\mu s$', '$t=2000\mu s$'});
txb = {'$t=38.9\mu s$', '$t=208.4\mu s$', '$t=377.9\mu s$', '$t=547.7\mu s$', '$t=716.3\mu s$',...
       '$t=886.2\mu s$', '$t=1057.7\mu s$', '$t=1223.5\mu s$', '$t=2000\mu s$'};
tvb = {'$t=0.1\mu s$', '$t=114\mu s$', '$t=283.1\mu s$', '$t=452.6\mu s$', '$t=620.6\mu s$',...
       '$t=791.6\mu s$', '$t=959.5\mu s$', '$t=1135.6\mu s$', '$t=2000\mu s$'};

% positions A
subplot(2,2,1)
hold on
x = -0.2:0.001:0.15;
ci = 1;
for i = 1:2:(nA-1)
  disp(i)
  plot(x, pdf(fitdist(xA(i,:)', 'Normal'), x), 'color', [0 0 1]./ci,...
    'DisplayName',char(txa(i)), 'linestyle', '-', 'linewidth', 3/ci)
  ci = ci + 1;
end
ci = 1;
for i = 2:2:(nA-1)
  disp(i)
  plot(x, pdf(fitdist(xA(i,:)', 'Normal'), x), 'color', [1 0 0]./ci,...
    'DisplayName',char(txa(i)), 'linestyle', '-', 'linewidth', 3/ci)
  ci = ci + 1;
end
plot(x, pdf(fitdist(xA(end,:)', 'Normal'), x), 'color', 'k',...
    'DisplayName',char(txa(end)), 'linestyle', ':', 'linewidth', 2)
set(gca, 'fontsize', 14)
title('Distribution of positions for \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'EastOutside', 'interpreter', 'latex')
ylim([0 27])
xlim([-0.13 0.17])

% velocities A
subplot(2,2,2)
hold on
x = -3:0.01:2.5;
ci = 1;
for i = 1:2:(nA-1)
  disp(i)
  plot(x, pdf(fitdist(vA(i,:)', 'Normal'), x), 'color', [0 0 1]./ci,...
    'DisplayName',char(tva(i)), 'linestyle', '-', 'linewidth', 3/ci)
  ci = ci + 1;
end
ci = 1;
for i = 2:2:(nA-1)
  disp(i)
  plot(x, pdf(fitdist(vA(i,:)', 'Normal'), x), 'color', [1 0 0]./ci,...
    'DisplayName',char(tva(i)), 'linestyle', '-', 'linewidth', 3/ci)
  ci = ci + 1;
end
plot(x, pdf(fitdist(vA(end,:)', 'Normal'), x), 'color', 'k',...
    'DisplayName',char(tva(end)), 'linestyle', ':', 'linewidth', 2)
set(gca, 'fontsize', 14)
title('Distribution of velocities for \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'EastOutside', 'interpreter', 'latex')
ylim([0 1.75])
xlim([-2.5 2.3])

% positions B
subplot(2,2,3)
hold on
x = -0.15:0.001:0.2;
ci = 1;
for i = 1:2:(nB-1)
  disp(i)
  plot(x, pdf(fitdist(xB(i,:)', 'Normal'), x), 'color', [0 0 1]./ci,...
    'DisplayName',char(txb(i)), 'linestyle', '-', 'linewidth', 4/ci)
  ci = ci + 1;
end
ci = 1;
for i = 2:2:(nB-1)
  disp(i)
  plot(x, pdf(fitdist(xB(i,:)', 'Normal'), x), 'color', [1 0 0]./ci,...
    'DisplayName',char(txb(i)), 'linestyle', '-', 'linewidth', 4/ci)
  ci = ci + 1;
end
plot(x, pdf(fitdist(xB(end,:)', 'Normal'), x), 'color', 'k',...
    'DisplayName',char(txb(end)), 'linestyle', ':', 'linewidth', 2)
set(gca, 'fontsize', 14)
title('Distribution of positions for \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'EastOutside', 'interpreter', 'latex')
ylim([0 27])
xlim([-0.13 0.17])

% velocities B
subplot(2,2,4)
hold on
x = -3:0.01:3;
ci = 1;
for i = 1:2:(nB-1)
  disp(i)
  plot(x, pdf(fitdist(vB(i,:)', 'Normal'), x), 'color', [0 0 1]./ci,...
    'DisplayName',char(tvb(i)), 'linestyle', '-', 'linewidth', 5-ci)
  ci = ci + 1;
end
ci = 1;
for i = 2:2:(nB-1)
  disp(i)
  plot(x, pdf(fitdist(vB(i,:)', 'Normal'), x), 'color', [1 0 0]./ci,...
    'DisplayName',char(tvb(i)), 'linestyle', '-', 'linewidth', 5-ci)
  ci = ci + 1;
end
plot(x, pdf(fitdist(vB(end,:)', 'Normal'), x), 'color', 'k',...
    'DisplayName',char(tvb(end)), 'linestyle', ':', 'linewidth', 2)
set(gca, 'fontsize', 14)
title('Distribution of velocities for \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'EastOutside', 'interpreter', 'latex')
xlim([-2.5 2.3])
ylim([0 1.75])

saveas(fig2, 'task3.png')
setPDFsize;
saveas(fig2, 'task3.pdf')

