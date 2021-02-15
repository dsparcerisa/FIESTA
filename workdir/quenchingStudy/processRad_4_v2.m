clear all

%% Common values
dosePctError = 0.03; % 3%
doseMinError = 0.2; % Gy
SSiP = 80; % pixels
radius_cm = 0.05; % For mean value
sigFactor = 2.5;
distanciaBase = 5.3; % cm
deltas = -0.5:0.001:0.5; % For the Micke algorithm

%% Values specific to this run
I_inicial = 44;
I_final = 43;
Ierror = 2.4; 
NValidPoints = 11;
logName = 'irrLog_2020_07_08_15_07_50.log';
logPath = '/Users/dani/Documents/FIESTA/logs/2020_07_08';
filmType = 'EBT3unl';
E0 = 4; % MeV
radName = 'Rad4';
layerThickness_um = 14; % EBT3unl - automatizar
substrateThickness_um = 0; % EBT3unl - automatizar

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

% Posición 7
warning('Select position 7');
[X7,Y7] = getMaxCoordinates(I);

% Posición 21
warning('Select position 21');
[X21,Y21] = getMaxCoordinates(I);

% Posición 27
warning('Select position 27');
[X27,Y27] = getMaxCoordinates(I);


%% Do the processing
% XY vectors
Yvec = [X3-X1, Y3-Y1] / 4;
Xvec = [X21-X1, Y21-Y1] / 8 + [X27-X7, Y27-Y7] / 8; % Media de los dos a igual Z
nX = norm(Xvec);

% Calculo de thetaY
thetaY = atand(1-(norm(Yvec) / nX)) % el ángulo vertical es casi cero

% Calculo de thetaX mediante punto 7
dist17 = norm([X1-X7, Y1-Y7])
thetaX1 = atand((dist17/6/nX) - 1/3)

% Calculo de thetaX mediante punto 27
dist127 = norm([X1-X27, Y1-Y27])
thetaX2 = atand((dist127/6/nX) - 1)

thetaX = 0.5*(thetaX1 + thetaX2)

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
Zpositions_cm = [0 1 2 3 2 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10]'; % For incorrect Z2 in first round
% Zpositions_cm = [0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10]';
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

%% Si es satisfactorio, salvar
saveFile = [radName '.mat'];
save(saveFile);

%% Para comenzar desde aquí (específico de Rad1)
clearvars -except radName
saveFile = [radName '.mat']
load(saveFile);

%% SELECCIÓN DE SIGMAS
[finalSigmas, dfinalSigmas] = findSigmasInRC(allI,NValidPoints,CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits)

%% Plot
figure
Zpos = (0:10) + distanciaBase;
errorbar(Zpos, finalSigmas, dfinalSigmas, 'b.')
grid on

%% Corrección de intensidad
I_values = linspace(I_inicial, I_final, Npoints);
t1_s = 0.05;
Q_pC = I_values * t1_s;
NProt = Q_pC / 1.60217662e-7;
delta_NProt = Ierror ./ I_values .* NProt;
relI = I_values ./ I_values(1);

%% Process all elements in allI
radius_pixels = radius_cm*pixCM;
sigmaX_mm = nan(Npoints,1);
sigmaY_mm = nan(Npoints,1);
deltasigmaX_mm = nan(Npoints,1);
deltasigmaY_mm = nan(Npoints,1);
meanDoses_3mm = nan(Npoints,1);
stdDoses_3mm = nan(Npoints,1);
xcenters = zeros(Npoints,1);
ycenters = zeros(Npoints,1);
Zpos = (0:10) + distanciaBase;
Zpos = Zpos ./ cosd(thetaX);

%% Considerando doble irradiación
for i=1:NValidPoints  
    d1 = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    d2 = getDoseMicke(double(allI{i+11}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    d2.data = 0.2 * d2.data * I_values(i) / I_values(i+11);
    
    [~, xc1, yc1, sX1, sY1] = meanAndCenterMass(-d1.data,radius_pixels);
    [~, xc2, yc2, sX2, sY2] = meanAndCenterMass(-d2.data,radius_pixels);
   
    if i==3
        d3 = getDoseMicke(double(allI{5}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
        d3.data = d3.data * I_values(3) / I_values(5);
        [~, xc3, yc3, sX3, sY3] = meanAndCenterMass(-d3.data,radius_pixels);
        d13 = merge2Gaussians(d1.data, d3.data, [xc1 yc1], [xc3 yc3], (sX1+sY1)/2, (sX3+sY3)/2, sigFactor);
        allD{3} = merge2Gaussians(d13, d2.data, [0 0], [xc2 yc2], (sX1+sY1)/2, (sX2+sY2)/2, sigFactor);
    elseif i==5
        allD{5} = merge2Gaussians(zeros(size(d2.data)), d2.data, [0 0], [xc2 yc2], (sX2+sY2)/2, (sX2+sY2)/2, sigFactor);
    else
        allD{i} = merge2Gaussians(d1.data, d2.data, [xc1 yc1], [xc2 yc2], (sX1+sY1)/2, (sX2+sY2)/2, sigFactor);
    end
    % pause
end
%figure
%plotNimages(allD)

[~, XvalC, YvalC] = merge2Gaussians(d1.data, d2.data, [xc1 yc1], [xc2 yc2], (sX1+sY1)/2, (sX2+sY2)/2, sigFactor);
XvalC_mm = 10 * XvalC / pixCM;
YvalC_mm = 10 * YvalC / pixCM;
distance_mm = sqrt(XvalC_mm.^2 + YvalC_mm.^2);
maskRadius = distance_mm <= (radius_cm*10);

%% 
for i=1:NValidPoints 
    centeredDose = allD{i};
    meanDoses_3mm(i) = mean(centeredDose(maskRadius));
    stdDoses_3mm(i) = std(centeredDose(maskRadius));
    try       
        [sigmaX, sigmaY, deltasigmaX, deltasigmaY] = getSigmas(XvalC_mm,YvalC_mm,centeredDose,dosePctError*centeredDose+doseMinError);
        sigmaX_mm(i) = sigmaX;
        sigmaY_mm(i) = sigmaY;
        deltasigmaX_mm(i) = deltasigmaX;
        deltasigmaY_mm(i) = deltasigmaY;
    catch
        sigmaX_mm(i) = nan;
        sigmaY_mm(i) = nan;
        deltasigmaX_mm(i) = nan;
        deltasigmaY_mm(i) = nan;
    end
end

%% Ver que todo está OK
figure(1);

sigmaX_pix = sigmaX_mm/10*pixCM;
sigmaY_pix = sigmaY_mm/10*pixCM;

for i = 1:NValidPoints
    subplot(4,3,i);
    imagesc(XvalC_mm, YvalC_mm, allD{i});
    title(sprintf('%i',i));
    hold on
    errorbar(xcenters(i), ycenters(i), sigmaY_mm(i), sigmaY_mm(i), sigmaX_mm(i), sigmaX_mm(i), 'bo');
end


%% Predicción del tamaño del spot
theMask=1:NValidPoints;
wX = (deltasigmaX_mm.^(-2));
wY = (deltasigmaY_mm.^(-2));
finalSigmasX = sigmaX_mm;
finalSigmasY = sigmaY_mm;
finalSigmas = (finalSigmasX.*wX + finalSigmasY.*wY)./(wX+wY);
dfinalSigmas = (wX+wY).^(-0.5);
figure(3)
errorbar(Zpos(theMask), finalSigmasX(theMask), deltasigmaX_mm(theMask), 'b.');
hold on
errorbar(Zpos(theMask), finalSigmasY(theMask), deltasigmaY_mm(theMask), 'r.');
errorbar(Zpos(theMask), finalSigmas(theMask), dfinalSigmas(theMask), 'k.-');
ylabel('Sigma (mm)')
xlabel('Distance (cm)');
title('XY averaged spot size');
grid on


%% Para calcular la eficiencia relativa mejor juntar todos los valores de todas las medidas
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
E0kap_vector = energyStoppingPowerKapton(E0, [0 kaptonThickness_cm]);
E0kap = E0kap_vector(2);

% TODO!: mover esto para que sea susceptible de filmType
airVecPos = 0:0.01:20;
[Eair_vec, ~, Sw_vec] = energyStoppingPower(E0kap, airVecPos);
validMask = ~isnan(Eair_vec);
Epos_preFilm = interp1(airVecPos(validMask), Eair_vec(validMask), Zpos);
Sw_preFilm = interp1(airVecPos(validMask), Sw_vec(validMask), Zpos);
Sw_rel = Sw_preFilm ./ Sw_preFilm(1);

Sw_preFilm = Sw_preFilm(theMask);
Epos_preFilm = Epos_preFilm(theMask);

%% Estudio de las dosis relativas esperadas según la energía:
fA = @(r,sigma) 1-exp(-r.^2./(2.*sigma.^2));
dfA = @(r,sigma, dsigma) (1-fA(r,sigma))*r.^2 ./ sigma.^3 .* dsigma;
intFractions = fA(radius_cm*10, finalSigmas');
dIntFractions = dfA(radius_cm*10, finalSigmas', dfinalSigmas');
Edep_inR_MeV = NProt(theMask).*Sw_preFilm(theMask).*intFractions(theMask).*layerThickness_um/10000;
dEdep_inR_MeV = Edep_inR_MeV .* sqrt( (dIntFractions(theMask) ./ intFractions(theMask)).^2 + ...
                (Ierror ./ I_values(theMask)).^2); % Faltaría el error en LET
cylRho = 1.0; % g/cm3
cylinderMass_kg = 0.001 * cylRho * (layerThickness_um/10000)*pi*radius_cm*radius_cm;
predictedDose_Gy = Edep_inR_MeV * 1.60218e-13 / cylinderMass_kg;
dpredictedDose_Gy = dEdep_inR_MeV * 1.60218e-13 / cylinderMass_kg;

measDose_Gy = meanDoses_3mm(theMask)';
dmeasDose_Gy = stdDoses_3mm(theMask)';
ratios = measDose_Gy ./ predictedDose_Gy;
dratios = ratios .* sqrt( (dpredictedDose_Gy ./ predictedDose_Gy).^2 + (dmeasDose_Gy / measDose_Gy).^2 );

%% Plottear resultados
figure(2)
subplot(3,1,1);
errorbar(Zpos(theMask), finalSigmasX(theMask), deltasigmaX_mm(theMask), 'b.');
hold on
errorbar(Zpos(theMask), finalSigmasY(theMask), deltasigmaY_mm(theMask), 'r.');
errorbar(Zpos(theMask), finalSigmas(theMask), dfinalSigmas(theMask), 'k.-');
legend({'X', 'Y', 'Merged'},'Location','SouthEast');
ylabel('Sigma (mm)')
xlabel('Distance (cm)');
title('Spot size');
grid on

subplot(3,1,2);
errorbar(Zpos(theMask), measDose_Gy, dmeasDose_Gy,'r.');  hold on
errorbar(Zpos(theMask), predictedDose_Gy, dpredictedDose_Gy,'b.');  hold on
legend('Measured', 'Predicted');
ylabel('Mean dose (Gy)')
xlabel('Distance (cm)');
grid on
title('Central dose');

subplot(3,1,3);
errorbar(Zpos(theMask), ratios, dratios,'k.');  hold on
ylabel('Ratio')
xlabel('Distance (cm)');
grid on
title('Ratio measured / predicted');
ylim([0 1.2]);

%% Plottear resultados en función de LET
figure(3)
errorbar(Sw_preFilm(theMask), ratios, dratios, 'k.');
grid on
xlabel('LET (MeV/cm)')
ylabel('RE');
ylim([0 1.2]);

