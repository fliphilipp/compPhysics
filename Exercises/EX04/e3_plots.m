%% Task 2
clc; clf; clear all;

data_Ax = importdata('trajAx.data');
data_Av = importdata('trajAv.data');

data_Bx = importdata('trajBx.data');
data_Bv = importdata('trajBv.data');

% RELAXATION TIME A
figure(1);
subplot(2,5,1)
plot(data_Ax(:,1))

subplot(2,5,2)
plot(data_Ax(:,2))

subplot(2,5,3)
plot(data_Ax(:,3))

subplot(2,5,4)
plot(data_Ax(:,4))

subplot(2,5,5)
plot(data_Ax(:,5))

subplot(2,5,6)
plot(data_Av(:,1))

subplot(2,5,7)
plot(data_Av(:,2))

subplot(2,5,8)
plot(data_Av(:,3))

subplot(2,5,9)
plot(data_Av(:,4))

subplot(2,5,10)
plot(data_Av(:,5))

% RELAXATION TIME B
figure(2);
subplot(2,5,1)
plot(data_Bx(:,1))

subplot(2,5,2)
plot(data_Bx(:,2))

subplot(2,5,3)
plot(data_Bx(:,3))

subplot(2,5,4)
plot(data_Bx(:,4))

subplot(2,5,5)
plot(data_Bx(:,5))

subplot(2,5,6)
plot(data_Bv(:,1))

subplot(2,5,7)
plot(data_Bv(:,2))

subplot(2,5,8)
plot(data_Bv(:,3))

subplot(2,5,9)
plot(data_Bv(:,4))

subplot(2,5,10)
plot(data_Bv(:,5))