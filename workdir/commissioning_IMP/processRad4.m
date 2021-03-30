clear all; close all
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
addpath '/Users/dani/Documents/FIESTA/workdir/quenchingStudy'
basePath = '/Users/dani/Google Drive/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/21_03_01 Commissioning IMP y celulas';
fullPath = [basePath filesep 'R4.tif'];
I  = imread(fullPath);
I = uint16(I);

[pixCM, maxBits] = getImgMetaInfo(fullPath);

%% Apply general cropping
warning('Apply general cropping');

I = imcrop(I);
close all
imshow(I);

%% Base processing data
deltas = -0.005:0.001:0.005;
load('/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_HF_EBT3.mat')
radius_max_cm = 0.03;
radius_max_pixels = radius_max_cm*pixCM;
radius_pocillo_cm = 0.2;
radius_pocillo_pixels = radius_pocillo_cm*pixCM;

%% Number of items

Npoints = input('Select number of subimages: ')
allI = {};

close all;

% Process all elements in allI
sigmaX_mm = nan(Npoints,1);
sigmaY_mm = nan(Npoints,1);
deltasigmaX_mm = nan(Npoints,1);
deltasigmaY_mm = nan(Npoints,1);
meanDoses_max = nan(Npoints,1);
stdDoses_max = nan(Npoints,1);
meanDoses_4mm = nan(Npoints,1);
stdDoses_4mm = nan(Npoints,1);
xcenters = nan(Npoints,1);
ycenters = nan(Npoints,1);
FWHMx_mm = nan(Npoints,1);
FWHMy_mm = nan(Npoints,1);

for i=1:Npoints
    warning('Select subimage %i / %i', i, Npoints);
    figure(1);
    allI{i} = imcrop(I);
    close(1)
    
    figure(2)
    [dose, varMat, dr, dg, db] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    
    % coger el verde que no satura
    %if (max(dose.data(:))>8)
    dose.data = 0.5*(dg + dr);
    %end
    
    try
        [~, ~, ~, ~, ~, meanvalueMax, ~, ~, ~, ~, stdDoseMax] = meanAndCenterMass(-dose.data,radius_max_pixels);       
        [~, ~, ~, ~, ~, meanvalue, ~, ~, ~, ~, stdDose] = meanAndCenterMass(-dose.data,radius_pocillo_pixels);
        [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter, FWHMx, FWHMy] = getSigmas2(10*dose.getAxisValues('X'),10*dose.getAxisValues('Y'),dose.data,0.05*dose.data+0.01);  
    catch
        warning('Could not process meanAndCenterMass, using nans');
        mask = nan;
        xcentre = nan;
        ycentre  = nan;
        sigmaX  = nan;
        sigmaY = nan;
        meanvalue = nan;
        stdDose = nan;
        meanvalueMax = nan;
        stdDoseMax = nan;        
        deltasigmaX = nan;
        deltasigmaY = nan;

    end    
    
    sigmaX_mm(i) = sigmaX;
    sigmaY_mm(i) = sigmaY;
    deltasigmaX_mm(i) = deltasigmaX;
    deltasigmaY_mm(i) = deltasigmaY;
    meanDoses_max(i) = meanvalueMax;
    stdDoses_max(i) = stdDoseMax;    
    meanDoses_4mm(i) = meanvalue;
    stdDoses_4mm(i) = stdDose;
    xcenters(i) = xCenter; %xcentre;
    ycenters(i) = yCenter; %ycentre;    
    FWHMx_mm(i) = FWHMx;
    FWHMy_mm(i) = FWHMy;    
    
    fprintf('SigmaX (mm): %3.3f +- %3.3f\n', sigmaX, deltasigmaX);
    fprintf('SigmaY (mm): %3.3f +- %3.3f\n', sigmaY, deltasigmaY);
    fprintf('FWHMx (mm): %3.3f (%3.2f sigmaX)\n', FWHMx, FWHMx/sigmaX);
    fprintf('FWHMy (mm): %3.3f (%3.2f sigmaY)\n', FWHMy, FWHMy/sigmaY);
    fprintf('Mean dose (top): %3.3f +- %3.3f\n', meanvalueMax, stdDoseMax);
    fprintf('Mean dose (r2mm): %3.3f +- %3.3f\n\n', meanvalue, stdDose);
    
    figure(3);
    imshow(allI{i});
    title(i)
    
    figure(4);
    subplot(1,2,1);
    dose.plotSlice;
    colorbar
    subplot(1,2,2);
    varMat.plotSlice;
    
    pause
end

finalSigmas = 0.5*(sigmaX_mm + sigmaY_mm);

%% Plot
figure
subplot(1,2,1);
dose.plotSlice
set(gca,'ColorScale','log')
caxis([1 13])
colorbar
subplot(1,2,2);
[C,h] = contour(dose.getAxisValues('X'), dose.getAxisValues('Y'), dg', [0.5 1 2 3 5 10 12]);
grid on
xlabel('X'); ylabel('Y');
caxis([1 13])
set(gca,'ColorScale','log')
clabel(C,h);
