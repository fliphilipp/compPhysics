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

%% plot
fig1 = figure(1); set(fig1, 'Position', [100, 10, 1200, 800]);

% RELAXATION TIME A
subplot(2,2,1)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Ax(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_A(:,1) + data_A(:,2), 'b:','linewidth',2);
plot(data_time, data_A(:,1) - data_A(:,2), 'b:','linewidth',2)
h2 = plot(data_time, data_A(:,1),'r-','linewidth', 3);
set(gca, 'fontsize', 14)
title('Positions for relaxation time \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend([h1 h2 h3], {'sample trajectories', ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation'});
set(lgd,'fontsize',14)

subplot(2,2,2)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Av(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_A(:,3) + data_A(:,4), 'b:','linewidth',2);
plot(data_time, data_A(:,3) - data_A(:,4), 'b:','linewidth',2)
h2 = plot(data_time, data_A(:,3),'r-','linewidth', 3);
set(gca, 'fontsize', 14)
title('Velocities for relaxation time \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend([h1 h2 h3], {'sample trajectories', ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation'});
set(lgd,'fontsize',14)

% RELAXATION TIME B
subplot(2,2,3)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Bx(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_B(:,1) + data_B(:,2), 'b:','linewidth',2);
plot(data_time, data_B(:,1) - data_B(:,2), 'b:','linewidth',2)
h2 = plot(data_time, data_B(:,1),'r-','linewidth', 3);
set(gca, 'fontsize', 14)
title('Positions for relaxation time \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend([h1 h2 h3], {'sample trajectories', ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation'});
set(lgd,'fontsize',14)

subplot(2,2,4)
hold on
for traj = 1:ntraj
  h1 = plot(data_time, data_Bv(:,traj),'Color', [0.7 0.7 0.7]);
end
h3 = plot(data_time, data_B(:,3) + data_B(:,4), 'b:','linewidth',2);
plot(data_time, data_B(:,3) - data_B(:,4), 'b:','linewidth',2)
h2 = plot(data_time, data_B(:,3),'r-','linewidth', 3);
set(gca, 'fontsize', 14)
title('Velocities for relaxation time \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend([h1 h2 h3], {'sample trajectories', ['mean (' num2str(nmean) ' trajectories)'], 'standard deviation'});
set(lgd,'fontsize',14)

saveas(fig1, 'task2.png')
setPDFsize;
saveas(fig1, 'task2.pdf')

%% Task 3 data import

clear all;
close all;

xA = importdata('xdistA.data');
vA = importdata('vdistA.data') .* 1000;
xB = importdata('xdistB.data');
vB = importdata('vdistB.data') .* 1000;
n = size(xA,1);

%% Task 3 plotting
fig1 = figure(1); set(fig1, 'Position', [100, 10, 1200, 800]);
cols = get(gca,'colororder');

subplot(2,2,1)
hold on
x = -0.15:0.001:0.15;
for i = 2:n
  plot(x, pdf(fitdist(xA(i,:)', 'Normal'), x), 'color', cols(i,:),...
    'DisplayName',['$t = ' num2str(i-1) '/f_0$'], 'linestyle', '--')
end
set(gca, 'fontsize', 14)
title('Distribution of positions for \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'northwest', 'interpreter', 'latex')

subplot(2,2,2)
hold on
x = -3:0.01:3;
for i = 2:n
  plot(x, pdf(fitdist(vA(i,:)', 'Normal'), x), 'color', cols(i,:),...
    'DisplayName',['$t = ' num2str(i-1) '/f_0$'], 'linestyle', '--')
end
set(gca, 'fontsize', 14)
title('Distribution of velocities for \tau = 48.5 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'northwest', 'interpreter', 'latex')

subplot(2,2,3)
hold on
x = -0.15:0.001:0.15;
for i = 2:n
  plot(x, pdf(fitdist(xB(i,:)', 'Normal'), x), 'color', cols(i,:),...
    'DisplayName',['$t = ' num2str(i-1) '/f_0$'], 'linestyle', '--')
end
set(gca, 'fontsize', 14)
title('Distribution of positions for \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('x [\mum]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'northwest', 'interpreter', 'latex')

subplot(2,2,4)
hold on
x = -3:0.01:3;
for i = 2:n
  plot(x, pdf(fitdist(vB(i,:)', 'Normal'), x), 'color', cols(i,:),...
    'DisplayName',['$t = ' num2str(i-1) '/f_0$'], 'linestyle', '--')
end
set(gca, 'fontsize', 14)
title('Distribution of velocities for \tau = 147.3 \mus','fontsize',20)
xlabel('Time [\mus]','fontsize',20)
ylabel('v [\mum/ms]','fontsize',20)
lgd = legend('show');
set(lgd, 'fontsize', 14, 'location', 'northwest', 'interpreter', 'latex')

saveas(fig1, 'task3.png')
setPDFsize;
saveas(fig1, 'task3.pdf')
