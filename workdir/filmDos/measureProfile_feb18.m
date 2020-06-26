%% Measure calibration profile
clear all; clc; close all
FIESTAConfig
Fcfg.filmPath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_02_18 - Activacion huevos Camp 2';
Fcfg.filmCalFile = '/Users/dani/Documents/FIESTA/filmDosimetry/CoefFitCMAM2_outdated.mat'
load(Fcfg.filmCalFile, 'CoefR', 'CoefG', 'CoefB');

%% Define files to load
Nfiles = 2;
filePaths = {fullfile(Fcfg.filmPath,'Rad-1b0001.tif'),fullfile(Fcfg.filmPath,'Rad-1a0001.tif')};
nCrops = [3 3];

%% Read metadata
[pixelsXcm, maxInt] = getImgMetaInfo(filePaths{1});


%% ITERATE
N_ITERATIONS = 4;
doseCentralValue = nan(sum(nCrops), N_ITERATIONS);
stdCentralValue = nan(sum(nCrops), N_ITERATIONS);
integralDose = nan(sum(nCrops), N_ITERATIONS);

for ii=1:N_ITERATIONS
%% Read and crop all files
allI = loadNcropFiles(filePaths, nCrops);

%% Plot all images
%figure(1);
%plotNimages(allI);

%% Measure doses
deltas = 0.8:0.002:1.2;
allD = {};
for i=1:numel(allI)
    dose = getDoseFromRC(allI{i}, CoefR, CoefG, CoefB, pixelsXcm, deltas);
    allD{i} = dose.data;
    %meanDoses(i) = mean2(dose.data);
    %stdDoses(i) = std(dose.data(:));
end
% figure(2);
% plotNimages(allD);
% colorbar

%% 1. Medir sigmas y dosis central 
N = numel(allI);
sigmasX = nan(N, 1);
sigmasY = nan(N, 1);
deltasigmasX = nan(N, 1);
deltasigmasY = nan(N, 1);

radiusForIntegral_mm = 5;
radiusForIntegral_voxels = radiusForIntegral_mm * pixelsXcm / 10;
radiusInMM = 0.3;

radiusInVoxels = radiusInMM * pixelsXcm / 10;

for i=1:N
    [mask, xcentre, ycentre, sigmaX, sigmaY, meanValue, mX, mY, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(-allD{i}, radiusInVoxels);
    [m1,         ~,       ~,      ~,      ~,        m2] = meanAndCenterMass(-allD{i}, radiusForIntegral_voxels);
    integralDose(i, ii) = m2*numel(m1);   
    
    sigmasX(i) = sigmaX;
    sigmasY(i) = sigmaY;
    deltasigmasX(i) = deltasigmaX;
    deltasigmasY(i) = deltasigmaY;        
    doseCentralValue(i,ii) = meanValue;
    stdCentralValue(i,ii) = stdDose;
    
end
sigmasX_mm(:,ii) = sigmasX / pixelsXcm * 10;
sigmasY_mm(:,ii) = sigmasY / pixelsXcm * 10;

deltasigmasX_mm(:,ii) = deltasigmasX / pixelsXcm * 10;
deltasigmasY_mm(:,ii) = deltasigmasY / pixelsXcm * 10;

end

save('contours_feb18.mat');

%% Final values
clear all
load('contours_feb18.mat');

meanIntegralDose = mean(integralDose, 2, 'omitnan');
errorIntegralDose = std(integralDose, [], 2, 'omitnan');

meanSigmaX = mean(sigmasX_mm, 2, 'omitnan')
meanSigmaY = mean(sigmasY_mm, 2, 'omitnan')
meanDoses = mean(doseCentralValue, 2)

errorSigmaX = sqrt(mean(deltasigmasX_mm,2,'omitnan').^2 + (std(sigmasX_mm,[],2,'omitnan')/sqrt(N_ITERATIONS)).^2)
errorSigmaY = sqrt(mean(deltasigmasY_mm,2,'omitnan').^2 + (std(sigmasY_mm,[],2,'omitnan')/sqrt(N_ITERATIONS)).^2)
errorDoses = sqrt(mean(stdCentralValue,2).^2 + (std(doseCentralValue,[],2)/sqrt(N_ITERATIONS)).^2)

%% Define Z values and make plots
z0 = 5.0;
zValues = z0 + (0:1.5:7.5);
[zOrdered, izO] = sort(zValues);

figure(3)
subplot(1,2,1);
errorbar(zValues, meanSigmaX, errorSigmaX, 'bo'); hold on
errorbar(zValues, meanSigmaY, errorSigmaY, 'rx')
grid on
xlabel('Air distance (cm)');
ylabel('Beam sigma (mm)');
legend('SigmaX', 'SigmaY');
title('Sigma');
subplot(1,2,2)
errorbar(zValues, meanDoses, errorDoses,'bo'); hold on
ylabel('Central dose (Gy)');
xlabel('Air distance (cm)');
title('Dose');
grid on

%% Output values
meanSigmaX(izO)
errorSigmaX(izO)
meanSigmaY(izO)
errorSigmaY(izO)
meanDoses(izO)
errorDoses(izO)
aa = integralDose(izO)' * 1e-4;
aa / aa(1)
errorIntegralDose(izO) * 1e-4 / aa(1)

%% Make plot to fit
%finalSigmas = (sigmasX_mm + sigmasY_mm) / 2;
%subplot(1,2,2);
hold off
errorbar(zValues, meanSigmaX, errorSigmaX, 'bo'); 
hold on
errorbar(zValues, meanSigmaY, errorSigmaY, 'ro'); 

%F1x = fit(zValues', meanSigmaX, 'poly2', 'Lower', [0 0 0], 'Upper', [1 1 0])
F2x = fit(zValues', meanSigmaX, 'poly2', 'Weight', errorSigmaX.^(-2))
F2y = fit(zValues', meanSigmaY, 'poly2', 'Weight', errorSigmaY.^(-2))

%F3x = fit(zValues', meanSigmaX, 'poly1')

hold on
%plot(F1x, 'b')
plot(F2x, 'b')
plot(F2y, 'r')
%plot(F3x, 'm')
grid on
legend('SigmaX', 'SigmaY');
xlabel('Air depth (cm)');
ylabel('Sigma (mm)');
%MCVal = 10*(0.0016*10.^2 + 0.0048*10);
%plot(zValues, MCVal, 'k')

%%
% sigmasPP = 0.16 * 0.064*zValues;
% plot(zValues, sigmasPP, 'r-')
% legend('Sigmas sin pepperpot', 'Fit', 'Fit sigmas con pepperpot')
xlabel('Z [cm]'); ylabel('Sigma [mm]'); title('Sin difusor');
set(gca, 'FontSize', 14)
grid on
subplot(1,2,1);
fullImage=imread(filePaths{1});
if size(fullImage,3)>3
    fullImage = fullImage(:,:,1:3);
end
imshow(fullImage)

