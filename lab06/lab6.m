clear all
close all
clc

opengl software
% load the data
load heater.mat

% remove mean values
ue = ue-mean(ue);
ye = ye-mean(ye);
uv = uv-mean(uv);
yv = yv-mean(yv);

% define datasets
ze = [ye, ue];
zv = [yv, uv];

% ARX identification
for nk=1:5 % loop on I/O delay
    for na=1:5 % loop on the order na (nb=na is assumed)
        nb = na;
        model = arx(ze, [na, nb, nk]);
        figure, resid(ze, model, 'Corr') % check the residuals on the estimation dataset
    end
    pause(), close all
end