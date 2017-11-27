%% % Task 1 plots

clc; clf; clear all;
nbins = 25;

% load data file
distData = importdata('task1.data');

%plot
figure(1);
subplot(2,2,1)
histogram(distData(1:10,1), nbins,'Normalization', 'probability');
title('N = 10');
%xlim([0 1]);

subplot(2,2,2)
histogram(distData(1:100,2),nbins,'Normalization', 'probability');
title('N = 100');
%xlim([0 1]);

subplot(2,2,3)
histogram(distData(1:1000,3),nbins,'Normalization', 'probability');
title('N = 1000');
%xlim([0 1]);

subplot(2,2,4)
histogram(distData(1:10000,4),nbins,'Normalization', 'probability');
title('N = 10000');
xlim([0 1]);


%% Task 2 plots

clc; clf; clear all;
nbins = 25;

% energies
% load the data file
distData = importdata('task2.data');

%plot
figure(1);
subplot(2,2,1)
histogram(distData(1:10,1), nbins,'Normalization', 'probability');
title('N = 10');
xlim([0 1]);

subplot(2,2,2)
histogram(distData(1:100,2),nbins,'Normalization', 'probability');
title('N = 100');
xlim([0 1]);

subplot(2,2,3)
histogram(distData(1:1000,3),nbins,'Normalization', 'probability');
title('N = 1000');
xlim([0 1]);

subplot(2,2,4)
histogram(distData(1:10000,4),nbins,'Normalization', 'probability');
title('N = 10000');
xlim([0 1]);


%% Task 3

clc; clf; clear all;
nbins = 25;

data = importdata('task3.data');

% plot
figure(1);
subplot(3,3,1);
histogram(data(:,1), nbins, 'Normalization','probability');
title('Distribution of x');
xlabel('Probability');
ylabel('x');

subplot(3,3,2);
histogram(data(:,2), nbins, 'Normalization','probability');
title('Distribution of y');
xlabel('Probability');
ylabel('y');

subplot(3,3,3);
histogram(data(:,3), nbins, 'Normalization','probability');
title('Distribution of z');
xlabel('Probability');
ylabel('z');

subplot(3,3,4);
plot(1:10000,data(:,1));
title('Values of x over time');
xlabel('N');
ylabel('x');

subplot(3,3,5);
plot(1:10000,data(:,2));
title('Values of y over time');
xlabel('N');
ylabel('y');

subplot(3,3,6);
plot(1:10000,data(:,3));
title('Values of z over time');
xlabel('N');
ylabel('z');

subplot(3,3,7:9);
histogram(data(:,4), nbins, 'Normalization','probability');
title('Probability distribution of p_n');
xlabel('p_n');
ylabel('Probability');

 
%% Task 4a
clc; clf; clear all;

data = importdata('MC.txt');
numlags = 1000;
corr = autocorr(data, numlags);

data = importdata('task4a.data');

figure(1);
subplot(2,1,1)
plot(1:length(data), corr(1:length(data)));
title('Auto-correlation values calculated with autocorr in MatLab');
xlabel('Step');
ylabel('Auto-correlation');
hold on;
plot(15,corr(15),'r*');
hold off;

subplot(2,1,2)
plot(data);
title('Auto-correlation function evaluated in C');
xlabel('Step');
ylabel('Auto-correlation');
hold on;
plot(15,data(15),'r*');
hold off;

%% Task 4b
clc; clf; clear all;

figure(1);
data = importdata('task4b.data');
plot(data(6:end));
ylabel('Statistical inefficiency (s)');
xlabel('Block size (B)');
axis square;
% corr = autocorr(data, numlags);
% plot(0:numlags,corr);

