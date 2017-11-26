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

%% Task 4b
clc; clf; clear all;

data = importdata('MC.txt');
numlags = 1000;

% data = importdata('task4b.data');
% plot(1:length(data), data);
% 
% [corr,lags] = xcorr(data);
% 
% corr = gather(corr);
% 
% corr = autocorr(data, numlags);
% 
% subplot(1,1,1)
% plot(0:numlags,corr);
