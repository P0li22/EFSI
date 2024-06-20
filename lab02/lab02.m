clear all
close all
clc

load sensor.mat;

i1 = 1; i2 = 1;
while z(i1) < 0.013
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

% least squares std form definition
y = z_ret;
phi = [Vz_ret, ones(N, 1)];

% least squares solution
p = phi \ y
Kt = inv(p(1))
V0 = -p(2)/p(1)

% EUI_inf
A = pinv(phi);
eps = 5e-4;
p_min_eui = [0; 0];
for k = 1:N
    p_min_eui = p_min_eui + A(:,k).*(y(k)-eps*sign(A(:, k)));
end
p_min_eui
p_max_eui = 2*p-p_min_eui

% FPS_inf
M = [phi; -phi];

