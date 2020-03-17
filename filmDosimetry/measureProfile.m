%% Measure calibration profile
clear all; clc; close all
FIESTAConfig
load(Fcfg.filmCalFile, 'CoefR', 'CoefG', 'CoefB');


%% Define files to load
Nfiles = 3;
filePaths = {fullfile(Fcfg.filmPath,'Rad12a-Cal0001.tif'), fullfile(Fcfg.filmPath,'Rad12b-Cal0001.tif'), fullfile(Fcfg.filmPath,'Rad12c-Cal0001.tif')};
nCrops = [4 4 2];

%% Read metadata
[pixelsXcm, maxInt] = getImgMetaInfo(filePaths{1});

%% Read and crop all files
allI = loadNcropFiles(filePaths, nCrops)

%% Plot all images
figure(1);
plotNimages(allI);

%% Measure doses
deltas = 0.99:0.01:1.01;
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
zValues = 4.4 + (0:1.5:6);
lowDoseIndex = [4 5 9 2 6];
highDoseIndex = [3 7 10 4 8];
figure(3)
subplot(1,2,1);
plot(zValues, sigmasX_mm(lowDoseIndex), 'bo'); hold on
plot(zValues, sigmasX_mm(highDoseIndex), 'bx'); hold on
plot(zValues, sigmasY_mm(lowDoseIndex), 'ro'); hold on
plot(zValues, sigmasY_mm(highDoseIndex), 'rx')
grid on
xlabel('Air distance (cm)');
ylabel('Beam sigma (mm)');
legend('SigmaX low dose', 'SigmaX high dose', 'SigmaY low dose', 'SigmaY high dose');
title('Sigma');
subplot(1,2,2)
plot(zValues, doseCentralValue(lowDoseIndex), 'bo'); hold on
plot(zValues, doseCentralValue(highDoseIndex), 'rx'); hold on
ylabel('Central dose (Gy)');
xlabel('Air distance (cm)');
legend('Low dose', 'High dose');
title('Dose');
grid on

%% Make plot to fit
finalSigmas = (sigmasX_mm(highDoseIndex) + sigmasY_mm(highDoseIndex)) / 2;
subplot(1,2,2);
hold off
plot(zValues, finalSigmas, 'ko'); 
F1 = fit(zValues', finalSigmas, 'poly2', 'Lower', [0 0 0], 'Upper', [1 1 0])
%F2 = fit(zValues', finalSigmas, 'poly2')
%F3 = fit(zValues', finalSigmas, 'poly1')

hold on
plot(F1, 'b')
%plot(F2,'b')
%plot(F3,'m')

%MCVal = 10*(0.0016*10.^2 + 0.0048*10);
%plot(zValues, MCVal, 'k')

% sigmasPP = 0.16 * 0.064*zValues;
% plot(zValues, sigmasPP, 'r-')
% legend('Sigmas sin pepperpot', 'Fit', 'Fit sigmas con pepperpot')
xlabel('Z [cm]'); ylabel('Sigma [mm]'); title('Sin pepperpot');
set(gca, 'FontSize', 14)
grid on
subplot(3,2,1);
fullImage=imread(filePaths{1});
if size(fullImage,3)>3
    fullImage = fullImage(:,:,1:3);
end
imshow(fullImage)
subplot(3,2,3);
fullImage=imread(filePaths{2});
if size(fullImage,3)>3
    fullImage = fullImage(:,:,1:3);
end
imshow(fullImage)
subplot(3,2,5);
fullImage=imread(filePaths{3});
if size(fullImage,3)>3
    fullImage = fullImage(:,:,1:3);
end
imshow(fullImage)

