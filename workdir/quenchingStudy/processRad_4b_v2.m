clear all

%% Common values
dosePctError = 0.03; % 3%
doseMinError = 0.2; % Gy
radius_cm = 0.03; % For mean value
sigFactor = 2.5;
distanciaBase = 5.3; % cm
deltas = -0.5:0.001:0.5; % For the Micke algorithm

%% Values specific to this run
I_inicial = 44;
I_final = 44;
Ierror = 4; 
NValidPoints = 11
SSiP_inicial = 80; % pixels
SSiP_final = 90; % pixels
SSiPi = round(linspace(SSiP_inicial,SSiP_final,NValidPoints));
SSiPi(10) = SSiP_final;
SSiPi(11) = SSiP_final;
SSiPi = [SSiPi SSiPi];

logName = 'irrLog_2020_07_08_15_29_08.log';
logPath = '/Users/dani/Documents/FIESTA/logs/2020_07_08';
filmType = 'EBT3';
E0 = 4; % MeV
radName = 'Rad4b';

% Layer thicknesses
if strcmp(filmType,'EBT3unl')
    layerThickness_um = 14; 
    substrateThickness_um = 0;
elseif strcmp(filmType,'EBT3')
    layerThickness_um = 28;
    substrateThickness_um = 125; 
elseif strcmp(filmType,'EBT2')
    layerThickness_um = 30;
    substrateThickness_um = 80; % Oversimplification
end
%% Load
calFile = ['/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_' filmType '.mat'];
load(calFile);
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept';
fullPath = [basePath filesep filmType filesep sprintf('%i',E0) 'MeV'];
imgPath1 = [fullPath filesep radName '.tif'];
imgPath2 = [fullPath filesep radName ' (2).tif'];
imgPath3 = [fullPath filesep radName ' (3).tif'];

[pixCM, maxBits] = getImgMetaInfo(imgPath1);

I1 = imread(imgPath1);
I2 = imread(imgPath2);
I3 = imread(imgPath3);

% Merge all three scans of the image
I = (1/3) * (double(I1) + double(I2) + double(I3));
I = uint16(I);

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

% Posición 21
warning('Select position 21');
[X21,Y21] = getMaxCoordinates(I);


%% Do the processing
% XY vectors
% ESPECIFICO DE LA 4B
Yvec = [X3-X1, Y3-Y1] / 4;
Xvec = [X21-X1, Y21-Y1] / 4; %
nX = norm(Xvec);

% Calculo de thetaY
thetaY = atand(1-(norm(Yvec) / nX)) % el ángulo vertical es casi cero

% Calculo de thetaX mediante punto 7
thetaX = 2; % porque lo conocemos

% Prueba vectores
xy2 = [X1, Y1] + 2*Yvec;
xy7 = [X1, Y1] + 2*Xvec;
xy21 = [X1, Y1] + 4*Xvec;

figure(1);
imshow(I); hold on

% Resto
Xpositions_cm = 3 - [3 2 3 3 -1 -2 -1]';
Ypositions_cm = 2.5 - [2.5 1.5 0.5 -1.5 2.5 1.5 0.5]' ;
Zpositions_cm = [0 1 2 2 0 1 2]'; 
Xshifts_cm = tand(thetaX)*Zpositions_cm;
PxPositions =  (Xpositions_cm+Xshifts_cm) * Xvec + Ypositions_cm * Yvec + [X1,Y1];

% Prueba
imshow(I); hold on
plot(PxPositions(:,1), PxPositions(:,2), 'bo');

%% Do the cropping and extension
Npoints = size(PxPositions, 1);
figure(2);
allI = {};

zeroDoseR = maxBits * (CoefR1(1) - CoefR1(2)/CoefR1(3));
zeroDoseG = maxBits * (CoefG1(1) - CoefG1(2)/CoefG1(3));
zeroDoseB = maxBits * (CoefB1(1) - CoefB1(2)/CoefB1(3));

for i = 1:Npoints
    SSiP = SSiPi(i);
    allI{i} = imcrop(I, [PxPositions(i,1)-SSiP/2 PxPositions(i,2)-SSiP/2 SSiP SSiP]);
    if any(size(allI{i},[1 2])<(SSiP+1))
        fullM = ones(SSiP+1,SSiP+1,3);
        fullM(:,:,1) = zeroDoseR;
        fullM(:,:,2) = zeroDoseG;
        fullM(:,:,3) = zeroDoseB;
        fullM(1:size(allI{i},1),1:size(allI{i},2),:) = allI{i};
        allI{i} = uint16(fullM);
    end
    subplot(4,6,i);
    imshow(allI{i});
    title(sprintf('%i',i));
end

%% Fijar allI
allI2 = {};
allI2{1} = allI{1};
allI2{2} = allI{2};
allI2{3} = allI{3};
allI2{5} = allI{4};
allI2{12} = allI{5};
allI2{13} = allI{6};
allI2{14} = allI{7};
allI2{22} = [];
allI = allI2;
clear allI2;
%% Si es satisfactorio, salvar
saveFile = [radName '.mat'];
save(saveFile);

%% Para comenzar desde aquí
clear all
NValidPoints = 3;
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
radName = 'Rad4b';
saveFile = [radName '.mat']
load(saveFile);

%% SELECCIÓN DE SIGMAS
[finalSigmas, dfinalSigmas] = findSigmasInRC(allI,NValidPoints,CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits)
figure(2)
Zpos = (0:10) + distanciaBase;
errorbar(Zpos(1:NValidPoints), finalSigmas, dfinalSigmas, 'b.')
grid on
%% Procesar usando la función
processAnyRad(radName)



