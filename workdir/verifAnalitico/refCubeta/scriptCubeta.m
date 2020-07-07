clear all; close all;

%% 1. Irradiar plan de calibración

planPath = 'plan_CONV_feb6_V2_cubeta.txt'; 
planStructure = readPlan(planPath);
% Corregir por intensidad medida realmente: (ver logbook)
planStructure.I = 0.090;

dz = 0.01; % cm
targetTh = 0.0014; % 14 um
targetSPR = 1; % water-equivalent

dxy = 0.01; % 0.1 mm
sizeX = 20;
sizeY = 20;
doseCanvas = createEmptyCG2D(dxy, sizeX, sizeY);
factorImuestra = 1;
airDepthAtPos0 = 4.6; %cm

% Datos del 6-feb (R14)
% F2x = 
% 
%      Linear model Poly2:
%      F2x(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =   -0.008123  (-0.01547, -0.0007705)
%        p2 =      0.3744  (0.2612, 0.4876)
%        p3 =     -0.1246  (-0.5373, 0.2882)
% 
% F2y = 
% 
%      Linear model Poly2:
%      F2y(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =   -0.007649  (-0.0195, 0.004198)
%        p2 =      0.3502  (0.1695, 0.5309)
%        p3 =    -0.06355  (-0.7179, 0.5909)

sigmaPolyX = [-0.008123 0.3744 -0.1246];
sigmaPolyY = [-0.007649 0.3502 -0.06355];

limits = [-2 1 -4.5 1]; 
dose = getDoseFromPlanV2(doseCanvas, planStructure, dz, targetTh, targetSPR, factorImuestra, airDepthAtPos0, sigmaPolyX, sigmaPolyY);
dose.crop(limits);
dose.data = flipud(dose.data);
dose.plotSlice; colorbar

%% 2. Cargar radiocrómica 

% Measure calibration profile
filmCalFile = '/Users/dani/Documents/FIESTA/filmDosimetry/CoefFitCMAM_HG.mat';
load(filmCalFile, 'CoefR', 'CoefG', 'CoefB');
filmPathRC = '19.tif';
[pixelsXcm, maxInt] = getImgMetaInfo(filmPathRC);
I = imread(filmPathRC);
I = I(:,:,(1:3)); % quick fix for imread error

% cropping
rct = [ 140.5100   11.5100  749.9800  370.9800];
I = imcrop(I, rct);
doseRC = getDoseFromRC(I, CoefR, CoefG, CoefB, pixelsXcm);
doseRC.plotSlice; colorbar
caxis([0 15]);

%% Ejemplo figura conjunta (mejorar mucho)
subplot(1,2,1);
dose.data = fliplr(dose.data);
dose.plotSlice; colorbar
caxis([0 15]);
subplot(1,2,2);
doseRC.plotSlice; colorbar
caxis([0 15]);
