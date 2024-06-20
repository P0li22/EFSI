clear all
close all
clc

load sensor.mat;

i1 = 1; i2 = 1;
while z(i1) < 0.01
    i1 = i1+1;
end

while z(i2) < 0.035
    i2 = i2+1;
end

% linear region
z_ret = z(i1:i2);
Vz_ret = Vz(i1:i2);
N = length(z_ret);
n = 2;

% plot OG data
figure(1)
plot(z, Vz, 'o')
hold on, grid on, zoom on
xline(0.01), xline(0.035)

% least squares std form definition
y = z_ret;
phi = [Vz_ret, ones(N, 1)];

% least squares solution
p = phi \ y
Kt = inv(p(1))
V0 = -p(2)/p(1)

% plot the solution
Vz0 = linspace(-10, 10, 1000);
z_hat = p(1)*Vz0 + p(2);
plot(z_hat, Vz0, 'r-')

% confidence intervals with known variance
% e = 5e-4;
% k = 2; % 95% confidence
% sigma_e = e/k; 
% var_e = sigma_e^2;
% var_p = var_e * ((phi' * phi) \ eye(size(phi, 2))); % var_e * inv(phi'*phi)
% sigma_p = sqrt(var_p);
% p_int = [p-k*diag(sigma_p), p+k*diag(sigma_p)];
% Kt_int = [1/p_int(1,2), 1/(p_int(1,1))]
% V0_int = -[p_int(2,2) / p_int(1,1), p_int(2,1) / p_int(1,2)]

%confidence interval with unknown variance
k = 2; % 95% confidence
var_e = (1/(N)) * (y'*(eye(N)-phi*inv(phi'*phi)*phi')*y);
sigma_e = sqrt(var_e);
var_p = var_e * ((phi' * phi) \ eye(size(phi, 2))); % var_e * inv(phi'*phi)
sigma_p = sqrt(var_p);
p_int = [p-k*diag(sigma_p), p+k*diag(sigma_p)];
Kt_int = [1/p_int(1,2), 1/(p_int(1,1))]
V0_int = -[p_int(2,2) / p_int(1,1), p_int(2,1) / p_int(1,2)]

% plot the interval
% z_hat = p(1)*Vz0 + p(2);
z_max = max(p_int(1, 1)*Vz0, p_int(1, 2)*Vz0) + p_int(2, 2);
z_min = min(p_int(1, 1)*Vz0, p_int(1, 2)*Vz0) + p_int(2, 1);
plot(z_max, Vz0, 'r--')
plot(z_min, Vz0, 'r--')

