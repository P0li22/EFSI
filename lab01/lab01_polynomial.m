clear all
close all
clc

load sensor.mat;
N = length(z);
n = 3; % polynomial degree

% plot OG data
figure(1)
plot(z, Vz, 'o')
hold on, grid on, zoom on
xline(0.01), xline(0.035)

p = polyfit(Vz, z, n)

% plot fitting polyonomial
Vz0 = linspace(-8, 8, 1000);
z_hat = p(1)*Vz0.^3+p(2)*Vz0.^2+p(3)*Vz0+p(4);
plot(z_hat, Vz0, 'r-')