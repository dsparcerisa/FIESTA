clear all; close all;

%% 1. Irradiar plan de calibración

planPath = 'plan_calibrationRC.txt'; 
planStructure = readPlan(planPath);
dz = 0.01; % cm
targetTh = 0.0014; % 14 um
targetSPR = 1; % water-equivalent

dxy = 0.1; % 0.1 mm
sizeX = 20;
sizeY = 20;
doseCanvas = createEmptyCG2D(dxy, sizeX, sizeY);
factorImuestra = 1;
airDepthAtPos0 = 4.4; %cm

% Datos del 4-feb (R3)
% % F2x = 
% 
%      Linear model Poly2:
%      F2x(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =     0.00758  (0.0003605, 0.0148)
%        p2 =      0.2116  (0.09237, 0.3308)
%        p3 =     -0.2501  (-0.7062, 0.2059)
% 
% F2y = 
% 
%      Linear model Poly2:
%      F2y(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =    0.002342  (-0.007479, 0.01216)
%        p2 =      0.2993  (0.1414, 0.4571)
%        p3 =     -0.5561  (-1.145, 0.03249)

sigmaPolyX = [0.00758 0.2116 -0.2501];
sigmaPolyY = [0.00234 0.2993 -0.5561];
limits = [-11 1 -4.5 -2.5]; 
dose = getDoseFromPlanV2(doseCanvas, planStructure, dz, targetTh, targetSPR, factorImuestra, airDepthAtPos0, sigmaPolyX, sigmaPolyY);
dose.crop(limits);
dose.plotSlice; colorbar

%% 2. Cargar radiocrómica 

% Measure calibration profile
filmCalFile = '/Users/dani/Documents/FIESTA/filmDosimetry/CoefFitCMAM_HG.mat';
load(filmCalFile, 'CoefR', 'CoefG', 'CoefB');
filmPathRC = 'cal0005.tif';
[pixelsXcm, maxInt] = getImgMetaInfo(filmPathRC);
I = imread(filmPathRC);
I = I(:,:,(1:3)); % quick fix for imread error

% cropping
rct = [97.51       28.51      288.98     1287.98];
I = imcrop(I, rct);
doseRC = getDoseFromRC(I, CoefR, CoefG, CoefB, pixelsXcm);
doseRC.plotSlice

%% "Truco" para simular el error
planStructure_error = planStructure;
planStructure_error.t_s(10) = planStructure_error.t_s(9);
dose_error = getDoseFromPlanV2(doseCanvas, planStructure_error, dz, targetTh, targetSPR, factorImuestra, airDepthAtPos0, sigmaPolyX, sigmaPolyY);
dose_error.crop(limits);

%% FlipLR
dose_error.data = flipud(dose_error.data);
%% Plot de ejemplo
subplot(2,1,1);
dose_error.plotSlice; colorbar
caxis([0 12]);
subplot(2,1,2);
doseRC.plotSlice; colorbar
caxis([0 12]);