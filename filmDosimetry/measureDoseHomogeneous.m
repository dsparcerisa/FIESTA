% Measure dose using calibration

clear all; clc; close all
FIESTAConfig
load(Fcfg.filmCalFile, 'CoefR', 'CoefG', 'CoefB');

%% Define files to load
Nfiles = 1;
filePaths = {fullfile(Fcfg.filmPath,'90001.tif')};
nCrops = 6;

%% Read metadata
[pixelsXcm, maxInt] = getImgMetaInfo(filePaths{1});

%% Read and crop all files
allI = loadNcropFiles(filePaths, nCrops)

%% Plot all files
plotNimages(allI);

%% Read dose
deltas = 0.99:0.001:1.01;
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

%% Plot
figure(1);
errorbar(meanDoses, stdDoses, 'r.');
ylabel('Measured dose (Gy)');
xlabel('Sample #');
grid on
set(gca, 'FontSize', 14)

%% TODO: Plot dose VS delivered charge / number of shots and normalize. 