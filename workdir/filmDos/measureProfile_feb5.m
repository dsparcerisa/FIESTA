%% Measure calibration profile
clear all; clc; close all
FIESTAConfig
Fcfg.filmPath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_02_04 - FLASH Camp 2';
Fcfg.filmCalFile = '/Users/dani/Documents/FIESTA/filmDosimetry/CoefFitCMAM_HG.mat'
load(Fcfg.filmCalFile, 'CoefR', 'CoefG', 'CoefB');

%% Define files to load
Nfiles = 3;
% Usamos el punto a z = z0 + 4.5
filePaths = {fullfile(Fcfg.filmPath,'7a0001.tif'),fullfile(Fcfg.filmPath,'7b0001.tif'),fullfile(Fcfg.filmPath,'7c0001.tif')};
nCrops = [2 2 2];

%% Read metadata
[pixelsXcm, maxInt] = getImgMetaInfo(filePaths{1});

disp('Contour LOW DOSE profiles');
% Read and crop all files
allI_LD = loadNcropFiles(filePaths, nCrops);

disp('Contour HIGH DOSE profiles');
% Read and crop all files
allI_HD = loadNcropFiles(filePaths, nCrops)

% Plot all images
figure(1);
plotNimages(allI_LD);
figure(10)
plotNimages(allI_HD);

deltas = 0.8:0.002:1.2;
allD_LD = {};
allD_HD = {};

for i=1:numel(allI_LD)
    dose = getDoseFromRC(allI_LD{i}, CoefR, CoefG, CoefB, pixelsXcm, deltas);
    allD_LD{i} = dose.data;
    
    dose = getDoseFromRC(allI_HD{i}, CoefR, CoefG, CoefB, pixelsXcm, deltas);
    allD_HD{i} = dose.data;    
    %meanDoses(i) = mean2(dose.data);
    %stdDoses(i) = std(dose.data(:));
end
figure(2);
plotNimages(allD_LD);
figure(20);
plotNimages(allD_HD);

save('contours_feb5.mat');

%% Buscar centros y sigmas --> Realinear y recropear automáticamente para que ambas sean idénticas

clear all
load('contours_feb5.mat');
fac = 7.5; % how many times HD is supposed to be larger than LD
HL_multiplier = 0.95:0.005:1.05;

allD_comb = {};
for i=1:numel(allI_LD)
    [~, xcentre_LD, ycentre_LD, sigmaX_LD, sigmaY_LD] = meanAndCenterMass(-allD_LD{i}, 100)
    [~, xcentre_HD, ycentre_HD, sigmaX_HD, sigmaY_HD] = meanAndCenterMass(-allD_HD{i}, 100)
    for j=1:numel(HL_multiplier)
        cropX_HL = HL_multiplier(j)*4*min([sigmaX_LD sigmaX_HD]) / 2;
        cropY_HL = HL_multiplier(j)*4*min([sigmaY_LD sigmaY_HD]) / 2;
        rect_LD = [ (0.95+rand()/10)*(xcentre_LD-cropX_HL) (0.95+rand()/10)*(ycentre_LD-cropY_HL) floor(2*cropX_HL)  floor(2*cropY_HL) ];
        rect_HD = [ (0.95+rand()/10)*(xcentre_HD-cropX_HL) (0.95+rand()/10)*(ycentre_HD-cropY_HL) floor(2*cropX_HL)  floor(2*cropY_HL) ];
        newLD = imcrop(allD_LD{i},rect_LD);
        newHD = imcrop(allD_HD{i},rect_HD);
        size(newLD)
        size(newHD)
        try
            allD_comb{i,j} = max(newLD, newHD ./ fac);
        catch
            warning('Skipped');
            i
            j
            allD_comb{i,j} = newLD;
        end
    end
end

%% 1. Medir sigmas y dosis central 
N_iter = numel(HL_multiplier);
N = numel(allI_LD);
sigmasX = nan(N, N_iter);
sigmasY = nan(N, N_iter);
deltasigmasX = nan(N, N_iter);
deltasigmasY = nan(N, N_iter);
doseCentralValue = nan(N, N_iter);
stdCentralValue = nan(N, N_iter);
radiusForIntegral_mm = 5;
radiusForIntegral_voxels = radiusForIntegral_mm * pixelsXcm / 10;
integralDose = nan(N, N_iter);
radiusInMM = 0.3;
radiusInVoxels = radiusInMM * pixelsXcm / 10;

for i=1:N
    for j=1:N_iter
    [mask, xcentre, ycentre, sigmaX, sigmaY, meanValue, mX, mY, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(-allD_comb{i,j}, radiusInVoxels);
    [m1,         ~,       ~,      ~,      ~,        m2] = meanAndCenterMass(-allD_comb{i,j}, radiusForIntegral_voxels);
    integralDose(i,j) = m2*numel(m1);
    sigmasX(i,j) = sigmaX;
    sigmasY(i,j) = sigmaY;
    deltasigmasX(i,j) = deltasigmaX;
    deltasigmasY(i,j) = deltasigmaY;        
    doseCentralValue(i,j) = meanValue;
    stdCentralValue(i,j) = stdDose;
    end
end
sigmasX_mm = sigmasX ./ pixelsXcm .* 10;
sigmasY_mm = sigmasY ./ pixelsXcm .* 10;

deltasigmasX_mm = deltasigmasX ./ pixelsXcm .* 10;
deltasigmasY_mm = deltasigmasY ./ pixelsXcm .* 10;


%% Take means
meanIntegralDose = mean(integralDose, 2, 'omitnan');
errorIntegralDose = std(integralDose, [], 2, 'omitnan');

meanSigmaX = mean(sigmasX_mm, 2, 'omitnan')
meanSigmaY = mean(sigmasY_mm, 2)
meanDoses = mean(doseCentralValue, 2)
meanSigmaY(6) = nan;
errorSigmaX = sqrt(mean(deltasigmasX_mm,2,'omitnan').^2 + (std(sigmasX_mm,[],2, 'omitnan')/sqrt(N_iter)).^2)
errorSigmaY = sqrt(mean(deltasigmasY_mm,2).^2 + (std(sigmasY_mm,[],2)/sqrt(N_iter)).^2)
errorDoses = sqrt(mean(stdCentralValue,2).^2 + (std(doseCentralValue,[],2)/sqrt(N_iter)).^2)


%% Define Z values and make plots
z0 = 4.9;
zValues = z0 + [0 4.5 1.5 6 3 7.5];
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
F2y = fit(zValues(1:5)', meanSigmaY(1:5), 'poly2', 'Weight', errorSigmaY(1:5).^(-2))

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

