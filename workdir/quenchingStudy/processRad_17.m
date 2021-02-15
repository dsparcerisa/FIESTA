clear all
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept';
fullPath = [basePath filesep 'EBT3' filesep '10MeV'];
radName = 'Rad17';
imgPath1 = [fullPath filesep radName '.tif'];
imgPath2 = [fullPath filesep radName ' (2).tif'];
imgPath3 = [fullPath filesep radName ' (3).tif'];

I1 = imread(imgPath1);
I2 = imread(imgPath2);
I3 = imread(imgPath3);

% Merge all three scans of the image

I = (1/3) * (double(I1) + double(I2) + double(I3));
I = uint16(I);


%% Log (pendiente cargar con código de Miguel)
logName = 'irrLog_2020_07_08_18_13_35.log';
logPath = '/Users/dani/Documents/FIESTA/logs/2020_07_08';
distanciaBase = 5.3; % cm

%% Info
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

%% Apply general cropping
warning('Apply general cropping');

I = imcrop(I);
close all
imshow(I);
%% Select positions
% Posición 1
warning('Select position 1');
[X1,Y1] = getMaxCoordinates(I);

% Posición 3
warning('Select position 3');
[X3,Y3] = getMaxCoordinates(I);

% Posición 27
warning('Select position 27');
[X27,Y27] = getMaxCoordinates(I);

%% Do the processing
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
SSiP = 80;

Npoints = size(PxPositions, 1);
figure(2);
allI = {};
for i = 1:Npoints
    allI{i} = imcrop(I, [PxPositions(i,1)-SSiP/2 PxPositions(i,2)-SSiP/2 SSiP SSiP]);
    subplot(4,6,i);
    imshow(allI{i});
    title(sprintf('%i',i));
end

%% Si es satisfactorio, salvar
saveFile = [radName '.mat'];
save(saveFile, 'allI');

%% Para comenzar desde aquí (específico de Rad1)
% clear all
% load('Rad3.mat');
% basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept';
% fullPath = [basePath filesep 'EBT3unl' filesep '4MeV'];
% radName = 'Rad2';
% imgPath1 = [fullPath filesep radName '.tif'];
% imgPath2 = [fullPath filesep radName ' (2).tif'];
% imgPath3 = [fullPath filesep radName ' (3).tif'];
% Npoints = 22;
% distanciaBase = 5.3;

%% Process all elements in allI
[pixCM, maxBits] = getImgMetaInfo(imgPath1);
deltas = -0.5:0.001:0.5;
load('/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_EBT3.mat')

radius_cm = 0.03;
radius_pixels = radius_cm*pixCM;

sigmaX_mm = nan(Npoints,1);
sigmaY_mm = nan(Npoints,1);
deltasigmaX_mm = nan(Npoints,1);
deltasigmaY_mm = nan(Npoints,1);
meanDoses_3mm = nan(Npoints,1);
stdDoses_3mm = nan(Npoints,1);
xcenters = nan(Npoints,1);
ycenters = nan(Npoints,1);
Zpos = (0:10) + distanciaBase;

% Todo: considerar la doble irradiación
for i = 1:Npoints   
    [dose, varMat] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    try
        [mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, ~, ~, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(-dose.data,radius_pixels);
        [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter] = getSigmas(10*dose.getAxisValues('X'),10*dose.getAxisValues('Y'),dose.data,0.05*dose.data+0.01);  
    catch
        warning('Could not process meanAndCenterMass, using nans');
        mask = nan;
        xcentre = nan;
        ycentre  = nan;
        sigmaX  = nan;
        sigmaY = nan;
        meanvalue = nan;
        deltasigmaX = nan;
        deltasigmaY = nan;
        stdDose = nan;
    end
%     sigmaX_mm(i) = 10*sigmaX/pixCM;
%     sigmaY_mm(i) = 10*sigmaY/pixCM;
%     deltasigmaX_mm(i) = 10*deltasigmaX/pixCM;
%     deltasigmaY_mm(i) = 10*deltasigmaY/pixCM;

    sigmaX_mm(i) = sigmaX;
    sigmaY_mm(i) = sigmaY;
    deltasigmaX_mm(i) = deltasigmaX;
    deltasigmaY_mm(i) = deltasigmaY;
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
    subplot(4,6,i);
    imshow(allI{i});
    title(sprintf('%i',i));
    hold on
    errorbar(xcenters(i), ycenters(i), sigmaY_pix(i), sigmaY_pix(i), sigmaX_pix(i), sigmaX_pix(i), 'bo');
end


%% Predicción del tamaño de spot
theMask=1:11
finalSigmasX = 0.5*(sigmaX_mm(theMask) + sigmaX_mm(11 + theMask));
finalSigmasY = 0.5*(sigmaY_mm(theMask) + sigmaY_mm(11 + theMask));
finalSigmas = 0.5*(finalSigmasX + finalSigmasY);
figure(3)
plot(Zpos(theMask), finalSigmasX(theMask), 'bo');
hold on
plot(Zpos(theMask), finalSigmasY(theMask), 'rx');
plot(Zpos(theMask), finalSigmas(theMask), 'kd-');
ylabel('Sigma (mm)')
xlabel('Distance (cm)');
title('XY averaged spot size');
grid on

%% Corrección de intensidad
Npoints = 22; % para aqui
I_inicial = 32;
I_final = 30
Ierror = 2.6; 
I_values = linspace(I_inicial, I_final, Npoints+5);
I_values = I_values(6:end); % Nuevo para los planes con punto cero
relI = I_values ./ I_inicial
meanDoses_3mm_corr = relI' .* meanDoses_3mm;
stdDoses_3mm_corr = sqrt((relI' .* stdDoses_3mm).^2 + (Ierror/I_inicial .* meanDoses_3mm).^2);


%% Para calcular la eficiencia relativa mejor juntar todos los valores de todas las medidas
E0 = 10;
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
E0kap_vector = energyStoppingPowerKapton(E0, [0 kaptonThickness_cm]);
E0kap = E0kap_vector(2);

airVecPos = 0:0.01:20;
[Eair_vec, ~, Sw_vec] = energyStoppingPower(E0kap, airVecPos);
validMask = ~isnan(Eair_vec);
Epos_preFilm = interp1(airVecPos(validMask), Eair_vec(validMask), Zpos);
Sw_preFilm = interp1(airVecPos(validMask), Sw_vec(validMask), Zpos);
Sw_rel = Sw_preFilm ./ Sw_preFilm(1);

% Hacer pasar por 50um de agua -->
Epos_preActiveLayer = nan(size(Epos_preFilm));
Sw_preActiveLayer = nan(size(Epos_preFilm));
Sw_rel = Sw_preFilm ./ Sw_preFilm(1);

% Zpos = (0:10) + 5.3; % Repetimos aquí

coatingThickness_cm = 125*1e-4;
for i=1:numel(Epos_preFilm)
    [e_temp, Sw_temp] = energyStoppingPowerWater(Epos_preFilm(i), [0 coatingThickness_cm coatingThickness_cm]);
    Epos_preActiveLayer(i) = e_temp(2);
    Sw_preActiveLayer(i) = Sw_temp(3);
end

Epos_preActiveLayer'
Sw_preActiveLayer'

%% Estudio de las dosis relativas esperadas según la energía:
fA = @(r,sigma) 1-exp(-r.^2./(2.*sigma.^2));
intFractions = fA(radius_cm*10, finalSigmas);
intRelFractions = intFractions ./ intFractions(1);


%% Plottear resultados
figure(2)
subplot(2,1,1);
errorbar(Zpos(theMask), sigmaX_mm(theMask), deltasigmaX_mm(theMask),'r.'); hold on
errorbar(Zpos(theMask), sigmaY_mm(theMask), deltasigmaY_mm(theMask),'g.');
errorbar(Zpos(theMask), sigmaX_mm(11+theMask), deltasigmaX_mm(11+theMask),'b.');
errorbar(Zpos(theMask), sigmaY_mm(11+theMask), deltasigmaY_mm(11+theMask),'m.');
%legend('X, small dose', 'Y small dose');
legend({'X, small dose', 'Y small dose', 'X large dose', 'Y large dose'},'Location','SouthEast');
ylabel('Sigma (mm)')
xlabel('Distance (cm)');
title('Spot size');
grid on

subplot(2,1,2);
% errorbar(Zpos(theMask), meanDoses_3mm(theMask), stdDoses_3mm(theMask),'r.'); hold on 
% errorbar(Zpos2(theMask), 0.2*meanDoses_3mm(11+theMask), 0.2*stdDoses_3mm(11+theMask),'b.');
errorbar(Zpos(theMask), meanDoses_3mm_corr(theMask), stdDoses_3mm_corr(theMask),'r.');  hold on
errorbar(Zpos(theMask), 0.2*meanDoses_3mm_corr(11+theMask), 0.2*stdDoses_3mm_corr(11+theMask),'b.');
plot(Zpos(theMask), meanDoses_3mm_corr(1).*intRelFractions.*Sw_rel(theMask)', 'k:');
legend('Small dose', 'Large dose / 5');
ylabel('Mean dose (Gy)')
xlabel('Distance (cm)');
grid on
title('Central dose');

