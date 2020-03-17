% 22 enero
% Measure dose using calibration

clear all; clc; close all
FIESTAConfig
load(Fcfg.filmCalFile, 'CoefR', 'CoefG', 'CoefB');

%% Define files to load
Nfiles = 1;
filePaths = {fullfile(Fcfg.filmPath,'17.tif')};
nCrops = 8;

%% Read metadata
[pixelsXcm, maxInt] = getImgMetaInfo(filePaths{1});

%% Read and crop all files
allI = loadNcropFiles(filePaths, nCrops)

%% Plot all files
plotNimages(allI);

%% Read dose
deltas = 0.6:0.001:1.2;
allD = {};
for i=1:numel(allI)
    dose = getDoseFromRC(allI{i}, CoefR, CoefG, CoefB, pixelsXcm, deltas);
    allD{i} = dose.data
    allD{i} = imgaussfilt(allD{i}, 5); % Filter if necessary
end
figure(2);
plotNimages(allD);

%% Medir sigmas y dosis central (diametro 0.5 mm)

sigmasX = nan(numel(allI), 1);
sigmasY = nan(numel(allI), 1);
doseCentralValue = nan(numel(allI), 1);
doseWell = nan(numel(allI), 1);

radiusInMM = 0.25;
radiusInVoxels = radiusInMM * pixelsXcm / 10;
wellRadiusInMM = 1;
wellRadiusInVoxels = wellRadiusInMM * pixelsXcm / 10;

for i=1:numel(allI)
        try
        [~, xcentre, ycentre, sigmaX, sigmaY, meanValue, mX, mY] = meanAndCenterMass(-allD{i}, radiusInVoxels);
        sigmasX(i) = sigmaX;
        sigmasY(i) = sigmaY;   
        doseCentralValue(i) = meanValue;
        [~, ~, ~, ~, ~, meanWellValue, ~, ] = meanAndCenterMass(-allD{i}, wellRadiusInVoxels);
        doseWell(i) = meanWellValue;
        catch
            doseCentralValue(i) = mean2(allD{i});
        end
end

sigmasX_mm = sigmasX / pixelsXcm * 10;
sigmasY_mm = sigmasY / pixelsXcm * 10;

subplot(1,2,1);
plot(sigmasX_mm, 'ro');
hold on
plot(sigmasY_mm, 'bo');
legend('SigmaX','SigmaY', 'Location', 'Northwest');
grid on

subplot(1,2,2);
plot(doseCentralValue, 'ro');
hold on
plot(doseWell, 'bo');
legend('doseMax','doseWell', 'Location', 'Northwest');
grid on

%% TODO: Sacar la dosis por disparo.

