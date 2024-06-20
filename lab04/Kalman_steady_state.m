clear all
close all
clc
rng('default');

load data.mat;
N = 4000;
N0 = 100;

%system definition
A = [0.96 0.5 0.27 0.28; -0.125 0.96 -0.08 -0.07; 0 0 0.85 0.97; 0 0 0 0.99];
B = [1 -1 2 1]';
C = [0 2 0 0];
D = 0;
Bv1 = sqrt(15)*[0.5 0 0 1]';
V1 = Bv1*Bv1';
V2 = 2000;
V12 = 0;
n = 4;

% initial state
x(:, 1) = [30 40 -70 -10]';
P = 0.5*eye(n);

% system simulation
for t = 1:N
    v1(:, t) = mvnrnd(zeros(1, n), V1)';
    v2(t) = sqrt(V2)*randn;
    x(:, t+1) = A*x(:, t) + B*u(t) + v1(:, t);
    y(t) = C*x(:, t) + D*u(t) + v2(t);
end

% Kalman 1-step predictor and Kalman filter (steady state)
x_h(:, 1) = zeros(n, 1);
sys = ss(A, [B, eye(n)], C, [D, zeros(1,n)], 1);
[kalman_model, K, P, K0] = kalman(sys, V1, V2, 0);

for t = 1:N
   y_h(t) = C * x_h(:, t);
   e(t) = y(t) - y_h(t);
   x_f(:, t) = x_h(:, t) + K0*e(t);
   y_f(t) = C*x_f(:, t);
   x_h(:, t+1) = A*x_h(:, t) + B*u(t) + K*e(t);
end

for k = 1:4
        RMSE_x_h(k) = norm(x(k, N0+1:N) - x_h(k, N0+1:N))/sqrt(N - N0);
        RMSE_x_f(k) = norm(x(k, N0+1:N) - x_f(k, N0+1:N))/sqrt(N - N0);
end

RMSE_y_h = norm(y(N0+1:N) - y_h(N0+1:N))/sqrt(N - N0)
RMSE_y_f = norm(y(N0+1:N) - y_f(N0+1:N))/sqrt(N - N0)
RMSE_x_h, RMSE_x_f

T = 1:N;

for k=1:n
    figure, plot(T,x(k,1:N),'g', T,x_h(k,1:N),'r-.', T,x_f(k,1:N),'b--'),
    title(['State x_',num2str(k),'(t)']), legend('System S','Predictor K','Filter F')
end

figure, plot(T,y(1:N),'g', T,y_h(1:N),'r-.', T,y_f(1:N),'b--'),
title('Output y(t)'), legend('System S','Predictor K','Filter F')


 