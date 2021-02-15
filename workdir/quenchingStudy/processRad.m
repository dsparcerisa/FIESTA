clear all
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM';

%fullPath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept/EBT2/10MeV';
fullPath =  '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept/EBT3unl/3MeV'

%imgPath = [fullPath filesep 'Rad18.tif'];
imgPath = [fullPath filesep 'Rad2.tif'];
I = imread(imgPath);
I = imcrop(I);
% El plan es:
% # X (cm), Y (cm), Z(cm), t (s)
% 5, 0, 0, 0.05
% 5, 0, -2, 0.05
% 5, 0, -4, 0.05
% 5, 0, -6, 0.05
% 5, 0, -8, 0.05
% 3, 2.5, 0, 0.05
% 2, 1.5, -1, 0.05
% 3, 0.5, -2, 0.05
% 2, -0.5, -3, 0.05
% 3, -1.5, -4, 0.05
% 2, -2.5, -5, 0.05
% 1, 2.5, -6, 0.05
% 0, 1.5, -7, 0.05
% 1, 0.5, -8, 0.05
% 0, -0.5, -9, 0.05
% 1, -1.5, -10, 0.05
% -1, 2.5, 0, 0.25
% -2, 1.5, -1, 0.25
% -1, 0.5, -2, 0.25
% -2, -0.5, -3, 0.25
% -1, -1.5, -4, 0.25
% -2, -2.5, -5, 0.25
% -3, 2.5, -6, 0.25
% -4, 1.5, -7, 0.25
% -3, 0.5, -8, 0.25
% -4, -0.5, -9, 0.25
% -3, -1.5, -10, 0.25

% Posición 1
[X1,Y1] = getMaxCoordinates(I);

% Posición 3
[X3,Y3] = getMaxCoordinates(I);

% Posición 27
[X27,Y27] = getMaxCoordinates(I);

% XY vectors
Yvec = [X3-X1, Y3-Y1] / 4;
Xvec = [X27-X1, Y27-Y1] / 6;

% Prueba vectores
xy2 = [X1, Y1] + 2*Yvec;
xy7 = [X1, Y1] + 2*Xvec;
xy21 = [X1, Y1] + 4*Xvec;

figure(1);
imshow(I); hold on
plot([X1 X3 X27 xy2(1) xy7(1) xy21(1)], [Y1 Y3 Y27 xy2(2) xy7(2) xy21(2)] , 'bo');

% Resto
Xpositions_cm = 3 - [3 2 3 2 3 2 1 0 1 0 1 -1 -2 -1 -2 -1 -2 -3 -4 -3 -4 -3]';
Ypositions_cm = 2.5 - [2.5 1.5 0.5 -0.5 -1.5 -2.5 2.5 1.5 0.5 -0.5 -1.5 2.5 1.5 0.5 -0.5 -1.5 -2.5 2.5 1.5 0.5 -0.5 -1.5]' ;

PxPositions =  Xpositions_cm * Xvec + Ypositions_cm * Yvec + [X1,Y1];

% Prueba
imshow(I); hold on
plot(PxPositions(:,1), PxPositions(:,2), 'bo');

%% Do the cropping
%Square size in pixels: 100 x 100;
SSiP = 100;

Npoints = size(PxPositions, 1);
figure(2);
allI = {};
for i = 1:Npoints
    allI{i} = imcrop(I, [PxPositions(i,1)-SSiP/2 PxPositions(i,2)-SSiP/2 SSiP SSiP]);
    subplot(4,6,i);
    imshow(allI{i});
    title(sprintf('%i',i));
end

%% Process all elements in allI
[pixCM, maxBits] = getImgMetaInfo(imgPath);
deltas = -0.5:0.001:0.5;
load('/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_HF_EBT3unl.mat')

radius_cm = 0.03;
radius_pixels = radius_cm*pixCM;

sigmaX_mm = nan(numel(I),1);
sigmaY_mm = nan(numel(I),1);
deltasigmaX_mm = nan(numel(I),1);
deltasigmaY_mm = nan(numel(I),1);
meanDoses_3mm = nan(numel(I),1);
stdDoses_3mm = nan(numel(I),1);
xcenters = nan(numel(I),1);
ycenters = nan(numel(I),1);
Zpos = 0:10;

% Todo: considerar la doble irradiación
for i = 1:Npoints

[dose, varMat] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
[mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, ~, ~, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(-dose.data,radius_pixels);
sigmaX_mm(i) = 10*sigmaX/pixCM;
sigmaY_mm(i) = 10*sigmaY/pixCM;
deltasigmaX_mm(i) = 10*deltasigmaX/pixCM;
deltasigmaY_mm(i) = 10*deltasigmaY/pixCM;
meanDoses_3mm(i) = meanvalue;
stdDoses_3mm(i) = stdDose;
xcenters(i) = xcentre;
ycenters(i) = ycentre;
end
%% Ver que todo está OK
figure(1);

sigmaX_pix = sigmaX_mm/10*pixCM;
sigmaY_pix = sigmaY_mm/10*pixCM;

for i = 1:Npoints
    allI{i} = imcrop(I, [PxPositions(i,1)-SSiP/2 PxPositions(i,2)-SSiP/2 SSiP SSiP]);
    subplot(4,6,i);
    imshow(allI{i});
    title(sprintf('%i',i));
    hold on
    errorbar(xcenters(i), ycenters(i), sigmaY_pix(i), sigmaY_pix(i), sigmaX_pix(i), sigmaX_pix(i), 'bo');
end

%% Plottear resultados
figure(2)
subplot(2,1,1);
errorbar(Zpos(1:Npoints), sigmaX_mm(1:Npoints), deltasigmaX_mm(1:Npoints),'r.'); hold on
errorbar(Zpos(1:Npoints), sigmaY_mm(1:Npoints), deltasigmaY_mm(1:Npoints),'g.');
%errorbar(Zpos, sigmaX_mm(12:22), deltasigmaX_mm(12:22),'b.');
%errorbar(Zpos, sigmaY_mm(12:22), deltasigmaY_mm(12:22),'m.');
legend('X, small dose', 'Y small dose');
%legend('X, small dose', 'Y small dose', 'X large dose', 'Y large dose');
ylabel('Sigma (mm)')
xlabel('Distance (cm)');
title('Spot size');

subplot(2,1,2);
errorbar(Zpos(1:Npoints), meanDoses_3mm(1:Npoints), stdDoses_3mm(1:Npoints),'r.'); hold on 
%errorbar(Zpos, meanDoses_3mm(12:22), stdDoses_3mm(12:22),'b.');
%legend('Small dose', 'Large dose');
ylabel('Mean dose (Gy)')
xlabel('Distance (cm)');
title('Central dose');