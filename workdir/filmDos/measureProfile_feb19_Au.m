%% Measure calibration profile
clear all; clc; close all
FIESTAConfig
Fcfg.filmPath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocrómicas/20_02_18 - Activacion huevos Camp 2';
Fcfg.filmCalFile = '/Users/dani/Documents/FIESTA/filmDosimetry/CoefFitCMAM2_outdated.mat'
load(Fcfg.filmCalFile, 'CoefR', 'CoefG', 'CoefB');


%% Define files to load
Nfiles = 1;
filePaths = {fullfile(Fcfg.filmPath,'Rad-Au-10001.tif')};
nCrops = 6;

%% Read metadata
[pixelsXcm, maxInt] = getImgMetaInfo(filePaths{1});

%% Read and crop all files
allI = loadNcropFiles(filePaths, nCrops)

%% Plot all images
figure(1);
plotNimages(allI);

%% Measure doses
deltas = 0.75:0.001:1.25;
allD = {};
meanDoses = nan(numel(allI), 1);
stdDoses = nan(numel(allI), 1);
for i=1:numel(allI)
    dose = getDoseFromRC(allI{i}, CoefR, CoefG, CoefB, pixelsXcm, deltas);
    allD{i} = dose.data;
    meanDoses(i) = mean2(dose.data);
    stdDoses(i) = std(dose.data(:));
end
figure(2);
plotNimages(allD);


%% 1. Medir sigmas y dosis central 
N = numel(allI);
sigmasX = nan(N, 1);
sigmasY = nan(N, 1);
doseCentralValue = nan(N, 1);

radiusInMM = 0.25;
radiusInVoxels = radiusInMM * pixelsXcm / 10;

for i=1:N
    [mask, xcentre, ycentre, sigmaX, sigmaY, meanValue, mX, mY] = meanAndCenterMass(-allD{i}, radiusInVoxels);
    sigmasX(i) = sigmaX;
    sigmasY(i) = sigmaY;
    doseCentralValue(i) = meanValue;
end
sigmasX_mm = sigmasX / pixelsXcm * 10;
sigmasY_mm = sigmasY / pixelsXcm * 10;

%% Define Z values and make plots
zValues = 4.4 + (0:1.5:7.5);
figure(3)
subplot(1,2,1);
plot(zValues, sigmasX_mm, 'bo'); hold on
plot(zValues, sigmasY_mm, 'rx')
grid on
xlabel('Air distance (cm)');
ylabel('Beam sigma (mm)');
legend('SigmaX', 'SigmaY');
title('Sigma');
subplot(1,2,2)
plot(zValues, doseCentralValue, 'bo'); hold on
ylabel('Central dose (Gy)');
xlabel('Air distance (cm)');
title('Dose');
grid on

%% Make plot to fit
finalSigmas = (sigmasX_mm(1:5) + sigmasY_mm(1:5)) / 2;
subplot(1,2,2);
hold off
plot(zValues(1:5), finalSigmas, 'ko'); 
%F1 = fit(zValues(1:5)', finalSigmas, 'poly2', 'Lower', [0 0 0], 'Upper', [1 1 0])
F2 = fit(zValues(1:5)', finalSigmas, 'poly2')
%F3 = fit(zValues', finalSigmas, 'poly1')

hold on
%plot(F1, 'b')
plot(F2,'b')
%plot(F3,'m')

%MCVal = 10*(0.0016*10.^2 + 0.0048*10);
%plot(zValues, MCVal, 'k')

% sigmasPP = 0.16 * 0.064*zValues;
% plot(zValues, sigmasPP, 'r-')
% legend('Sigmas sin pepperpot', 'Fit', 'Fit sigmas con pepperpot')
xlabel('Z [cm]'); ylabel('Sigma [mm]'); title('Con difusor de oro');
set(gca, 'FontSize', 14)
grid on
subplot(1,2,1);
fullImage=imread(filePaths{1});
if size(fullImage,3)>3
    fullImage = fullImage(:,:,1:3);
end
imshow(fullImage)

