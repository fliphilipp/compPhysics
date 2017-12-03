%% Task 2
clc; clf; clear all;

data_Ax = importdata('trajAx.data');
data_Av = importdata('trajAv.data');

data_Bx = importdata('trajBx.data');
data_Bv = importdata('trajBv.data');

% rows = length(data_Ax);
rows = 1000;

% RELAXATION TIME A
figure(1);
subplot(2,2,1)
plot(1:rows, data_Ax(1:rows,1),1:rows, data_Ax(1:rows,2),1:rows, data_Ax(1:rows,3),1:rows, data_Ax(1:rows,4),1:rows, data_Ax(1:rows,5))
title('Positions for relaxation time \tau = 48.5 \mus')
xlabel('Steps')
ylabel('x \mum')

subplot(2,2,2)
plot(1:rows, data_Av(1:rows,1),1:rows, data_Av(1:rows,2),1:rows, data_Av(1:rows,3),1:rows, data_Av(1:rows,4),1:rows, data_Av(1:rows,5))
title('Velocities for relaxation time \tau = 48.5 \mus')
xlabel('Steps')
ylabel('v \mum/ms')

% RELAXATION TIME B
subplot(2,2,3)
plot(1:rows, data_Bx(1:rows,1),1:rows, data_Bx(1:rows,2),1:rows, data_Bx(1:rows,3),1:rows, data_Bx(1:rows,4),1:rows, data_Bx(1:rows,5))
title('Positions for relaxation time \tau = 147.3 \mus')
xlabel('Steps')
ylabel('x \mum')

subplot(2,2,4)
plot(1:rows, data_Bv(1:rows,1),1:rows, data_Bv(1:rows,2),1:rows, data_Bv(1:rows,3),1:rows, data_Bv(1:rows,4),1:rows, data_Bv(1:rows,5))
title('Velocities for relaxation time \tau = 147.3 \mus')
xlabel('Steps')
ylabel('v \mum/ms')

% FOR AVERAGING
data_A = importdata('trajA_Average.data');
data_B = importdata('trajB_Average.data');

% rows = length(data_A);
% rows = 1000;

% RELAXATION TIME A
figure(2);
subplot(2,2,1)
std_p = data_A(1:rows,1) + data_A(1:rows,2);
std_m = data_A(1:rows,1) - data_A(1:rows,2);
plot(1:rows, data_A(1:rows,1), 1:rows, std_p,'r--', 1:rows, std_m,'r--')
title('Averaged positions for relaxation time \tau = 48.5 \mus')
xlabel('Steps')
ylabel('x \mum')

subplot(2,2,2)
std_p = data_A(1:rows,3) + data_A(1:rows,4);
std_m = data_A(1:rows,3) - data_A(1:rows,4);
plot(1:rows, data_A(1:rows,3), 1:rows, std_p,'r--', 1:rows, std_m,'r--')
title('Averaged velocities for relaxation time \tau = 48.5 \mus')
xlabel('Steps')
ylabel('v \mum/ms')

% RELAXATION TIME B
subplot(2,2,3)
std_p = data_B(1:rows,1) + data_B(1:rows,2);
std_m = data_B(1:rows,1) - data_B(1:rows,2);
plot(1:rows, data_B(1:rows,2), 1:rows, std_p,'r--', 1:rows, std_m,'r--')
title('Averaged positions for relaxation time \tau = 147.3 \mus')
xlabel('Steps')
ylabel('x \mum')

subplot(2,2,4)
std_p = data_B(1:rows,3) + data_B(1:rows,4);
std_m = data_B(1:rows,3) - data_B(1:rows,4);
plot(1:rows, data_B(1:rows,3), 1:rows, std_p,'r--', 1:rows, std_m,'r--')
title('Averaged velocities for relaxation time \tau = 147.3 \mus')
xlabel('Steps')
ylabel('v \mum/ms')
