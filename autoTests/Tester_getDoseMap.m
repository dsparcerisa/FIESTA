% Pruebas
% 1: Cuadrado de 1x1 mm en Z = 0 para E = 3 MeV y 1e7 protones
% 2: Lo mismo pero para Z = 10;
% 3: Lo mismo pero para E = 8 MeV;

% doseMap = getDoseMap(E0, z, dz, Nprot, targetTh, targetDens, targetSPR, N0)

clear all; close all

% Prueba 1
E0 = 3;
z = 10;
dz = 0.001;
dxy = 0.01;
Nprot = 6.25e7; % 1 pC
targetTh = 0.001; % 10 um = 10^-5 m = 10^-3 cm
targetSPR = 1;
sizeX = 11;
sizeY = 11;
N0 = createEmptyCG2D(dxy, sizeX, sizeY);
squareHalfSize = 0.005;
N0.data( (N0.getXindex(-squareHalfSize)):(N0.getXindex(squareHalfSize)), (N0.getYindex(-squareHalfSize)):(N0.getYindex(squareHalfSize))) = 1;
N0.data = N0.data ./ sum(N0.data(:));

figure('Position', [10 1100 500 400]);
N0.plotSlice
colorbar
title('N0');

tic; P1 = getDoseMap(E0, z, dz, Nprot, targetTh, targetSPR, N0); toc

figure('Position', [10 470 500 400]);
P1.plotSlice
colorbar
title('Dose');