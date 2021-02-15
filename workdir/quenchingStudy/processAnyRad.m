function [Sw_preFilm, dSw_preFilm, ratios, dratios] = processAnyRad(radName, radius_cm) %, I_inicial, I_final, Npoints, Ierror)
%% Load initial data
saveFile = [radName '.mat'];
load(saveFile, 'allI', 'E0', 'I_inicial', 'I_final', 'Ierror', 'NValidPoints', 'distanciaBase', 'filmType', ...
    'pixCM', 'maxBits', 'thetaX');
Npoints = numel(allI);
theMask = 1:NValidPoints;

%% Particularidades de la RC
if E0==3 || strcmp(radName, 'Rad5')
    hasCeroPt = 0;
else
    hasCeroPt = 1;
end

if E0==3 || strcmp(radName, 'Rad3') || strcmp(radName, 'Rad4') || strcmp(radName, 'Rad4b')
    error35 = 1;
else
    error35 = 0;
end

%% Calibration
calFile = ['/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_' filmType '.mat'];
load(calFile);
deltas = -0.5:0.001:0.5;

%% Process sigmas
sigmaFile = ['sigmas_' filmType '.mat'];
load(sigmaFile, 'sigmas');
switch(E0)
    case 3
        finalSigmas = sigmas(1:NValidPoints,1);
        dfinalSigmas = sigmas(1:NValidPoints,2);
    case 4
        finalSigmas = sigmas(1:NValidPoints,3);
        dfinalSigmas = sigmas(1:NValidPoints,4);
    case 5
        finalSigmas = sigmas(1:NValidPoints,5);
        dfinalSigmas = sigmas(1:NValidPoints,6);
    case 6
        finalSigmas = sigmas(1:NValidPoints,7);
        dfinalSigmas = sigmas(1:NValidPoints,8);
    case 8
        finalSigmas = sigmas(1:NValidPoints,9);
        dfinalSigmas = sigmas(1:NValidPoints,10);
    case 10
        finalSigmas = sigmas(1:NValidPoints,11);
        dfinalSigmas = sigmas(1:NValidPoints,12);
end

%% TEST
% finalSigmas = finalSigmas*0.8;
% warning('Artificially making sigmas 20% smaller');

%% Corrección de intensidad
if (hasCeroPt)
    I_values = linspace(I_inicial, I_final, Npoints+6); 
    I_values = I_values(7:end);
else
    I_values = linspace(I_inicial, I_final, Npoints);
end

t1_s = 0.05;
Q_pC = I_values * t1_s;
NProt = Q_pC / 1.60217662e-7;

%% Process all elements in allI
radius_pixels = radius_cm*pixCM;
meanDoses_3mm = nan(Npoints,1);
stdDoses_3mm = nan(Npoints,1);
xcenters = zeros(Npoints,1);
ycenters = zeros(Npoints,1);
Zpos = (0:10) + distanciaBase;
Zpos = Zpos ./ cosd(thetaX);
Zpos = Zpos(1:NValidPoints);
allD = {};

for i=1:NValidPoints
     d1 = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
     d2 = getDoseMicke(double(allI{i+11}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
     d2.data = 0.2 * d2.data * I_values(i) / I_values(i+11);   
     allD{i} = d1.data;
     allD{i+11} = d2.data;
     [~, xC, yC, ~, ~, meanvalue, ~, ~, ~, ~, stdDose] = meanAndCenterMass(-d1.data,radius_pixels);
     meanDoses_3mm(i) = meanvalue;
     stdDoses_3mm(i) = stdDose;
     xcenters(i) = xC;
     ycenters(i) = yC;
     [~, xC, yC, ~, ~, meanvalue, ~, ~, ~, ~, stdDose] = meanAndCenterMass(-d2.data,radius_pixels);
     meanDoses_3mm(i+11) = meanvalue;
     stdDoses_3mm(i+11) = stdDose;     
     xcenters(i+11) = xC;
     ycenters(i+11) = yC;     
end

finalMeanDoses = nan(NValidPoints,1);
dfinalMeanDoses = nan(NValidPoints,1);

for i = 1:NValidPoints
    if error35 && i==5
        finalMeanDoses(i) = meanDoses_3mm(i+11);
        dfinalMeanDoses(i) = stdDoses_3mm(i+11);
    elseif error35 && i==3
        [maxDose3, maxPos3] = max([meanDoses_3mm(3), meanDoses_3mm(5), meanDoses_3mm(14)]);
        allPos3 = [3 5 14];
        realPos3 = allPos3(maxPos3);
        finalMeanDoses(3) = meanDoses_3mm(realPos3);
        dfinalMeanDoses(3) = stdDoses_3mm(realPos3);               
    else
        if meanDoses_3mm(i) > meanDoses_3mm(i+11)
            finalMeanDoses(i) = meanDoses_3mm(i);
            dfinalMeanDoses(i) = stdDoses_3mm(i);
        else
            finalMeanDoses(i) = meanDoses_3mm(i+11);
            dfinalMeanDoses(i) = stdDoses_3mm(i+11);
        end
    end
   
end

% %% Considerando doble irradiación
% sigFactor = 4;
% for i=1:NValidPoints  
%     d1 = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
%     d2 = getDoseMicke(double(allI{i+11}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
%     d2.data = 0.2 * d2.data * I_values(i) / I_values(i+11);
%     
%     [~, xc1, yc1, sX1, sY1] = meanAndCenterMass(-d1.data,radius_pixels);
%     [~, xc2, yc2, sX2, sY2] = meanAndCenterMass(-d2.data,radius_pixels);
%    
%     % TODO: esto vale para la primera, corregir el resto
%     if i==3
%         d3 = getDoseMicke(double(allI{5}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
%         d3.data = d3.data * I_values(3) / I_values(5);
%         [~, xc3, yc3, sX3, sY3] = meanAndCenterMass(-d3.data,radius_pixels);
%         d13 = merge2Gaussians(d1.data, d3.data, [xc1 yc1], [xc3 yc3], (sX1+sY1)/2, (sX3+sY3)/2, sigFactor);
%         allD{3} = merge2Gaussians(d13, d2.data, [0 0], [xc2 yc2], (sX1+sY1)/2, (sX2+sY2)/2, sigFactor);
%     elseif i==5
%         allD{5} = merge2Gaussians(zeros(size(d2.data)), d2.data, [0 0], [xc2 yc2], (sX2+sY2)/2, (sX2+sY2)/2, sigFactor);
%     else
%         allD{i} = merge2Gaussians(d1.data, d2.data, [xc1 yc1], [xc2 yc2], (sX1+sY1)/2, (sX2+sY2)/2, sigFactor);
%     end
%     title(i)
% end
% %figure
% %plotNimages(allD)
% 
% % Hace un merge fake para hallar las coordenadas centradas y en mm.
% [~, XvalC, YvalC] = merge2Gaussians(d1.data, d2.data, [xc1 yc1], [xc2 yc2], (sX1+sY1)/2, (sX2+sY2)/2, sigFactor);
% XvalC_mm = 10 * XvalC / pixCM;
% YvalC_mm = 10 * YvalC / pixCM;
% distance_mm = sqrt(XvalC_mm.^2 + YvalC_mm.^2);
% maskRadius = distance_mm <= (radius_cm*10);
% 
% % Calcular valores medios
% for i=1:NValidPoints 
%     centeredDose = allD{i};
%     meanDoses_3mm(i) = mean(centeredDose(maskRadius));
%     stdDoses_3mm(i) = std(centeredDose(maskRadius));
% end

%% Ver que todo está OK
figure(1);
finalSigmas_pix = finalSigmas/10*pixCM;
for i = 1:NValidPoints
    subplot(4,3,i);
    imagesc(allD{i});
    title(sprintf('%i',i));
    hold on
    errorbar(xcenters(i), ycenters(i), finalSigmas_pix(i), finalSigmas_pix(i), finalSigmas_pix(i), finalSigmas_pix(i), 'bo');
end

%% Predicción del tamaño del spot
% theMask=1:NValidPoints;
% wX = (deltasigmaX_mm.^(-2));
% wY = (deltasigmaY_mm.^(-2));
% finalSigmasX = sigmaX_mm;
% finalSigmasY = sigmaY_mm;
% finalSigmas = (finalSigmasX.*wX + finalSigmasY.*wY)./(wX+wY);
% dfinalSigmas = (wX+wY).^(-0.5);
% figure(3)
% errorbar(Zpos(theMask), finalSigmasX(theMask), deltasigmaX_mm(theMask), 'b.');
% hold on
% errorbar(Zpos(theMask), finalSigmasY(theMask), deltasigmaY_mm(theMask), 'r.');
% errorbar(Zpos(theMask), finalSigmas(theMask), dfinalSigmas(theMask), 'k.-');
% ylabel('Sigma (mm)')
% xlabel('Distance (cm)');
% title('XY averaged spot size');
% grid on


%% Para calcular la eficiencia relativa mejor juntar todos los valores de todas las medidas
% kaptonThickness_um = 8;
% kaptonThickness_cm = 1e-4*kaptonThickness_um;
% E0kap_vector = energyStoppingPowerKapton(E0, [0 kaptonThickness_cm]);
% E0kap = E0kap_vector(2);
% 
% % TODO!: mover esto para que sea susceptible de filmType
% airVecPos = 0:0.01:20;
% [Eair_vec, ~, Sw_vec] = energyStoppingPower(E0kap, airVecPos);
% validMask = ~isnan(Eair_vec);
% Epos_preFilm = interp1(airVecPos(validMask), Eair_vec(validMask), Zpos);
% Sw_preFilm = interp1(airVecPos(validMask), Sw_vec(validMask), Zpos);
% Sw_rel = Sw_preFilm ./ Sw_preFilm(1);

% Load LET values
LETpath = '/Users/dani/Google Drive/UCM/01-Proyectos/25-TOPAS/LETScorer2/results';
LETFileName = sprintf('LETvalues_%s_%iMeV',filmType,E0);
fullLETPath = [LETpath filesep LETFileName];
load(fullLETPath);
theMask = 1:NValidPoints;
Sw_preFilm = 10*meanLET(theMask);
dSw_preFilm = 10*errLET(theMask);

%Sw_preFilm = Sw_preFilm(theMask);
% Epos_preFilm = Epos_preFilm(theMask);

%% Estudio de las dosis relativas esperadas según la energía:
if strcmp(filmType, 'EBT3unl')
    layerThickness_um = 14;
elseif strcmp(filmType, 'EBT3')
    layerThickness_um = 28;
elseif strcmp(filmType, 'EBT2')
    layerThickness_um = 30;
end

fA = @(r,sigma) 1-exp(-r.^2./(2.*sigma.^2));
dfA = @(r,sigma, dsigma) (1-fA(r,sigma))*r.^2 ./ sigma.^3 .* dsigma;
intFractions = fA(radius_cm*10, finalSigmas');
dIntFractions = dfA(radius_cm*10, finalSigmas', dfinalSigmas');
Edep_inR_MeV = NProt(theMask).*Sw_preFilm(theMask).*intFractions(theMask).*layerThickness_um/10000;
dEdep_inR_MeV = Edep_inR_MeV .* sqrt( (dIntFractions(theMask) ./ intFractions(theMask)).^2 + ...
                (Ierror ./ I_values(theMask)).^2 + (dSw_preFilm ./ Sw_preFilm).^2); 
cylRho = 1.2; % g/cm3
cylinderMass_kg = 0.001 * cylRho * (layerThickness_um/10000)*pi*radius_cm*radius_cm;
predictedDose_Gy = Edep_inR_MeV * 1.60218e-13 / cylinderMass_kg;
dpredictedDose_Gy = dEdep_inR_MeV * 1.60218e-13 / cylinderMass_kg;

measDose_Gy = finalMeanDoses';
dmeasDose_Gy = dfinalMeanDoses';
ratios = measDose_Gy ./ predictedDose_Gy;
dratios = ratios .* sqrt( (dpredictedDose_Gy ./ predictedDose_Gy).^2 + (dmeasDose_Gy / measDose_Gy).^2 );

%% Plottear resultados
figure(2)
subplot(3,1,1);
errorbar(Zpos, finalSigmas, dfinalSigmas, 'k.-');
ylabel('Sigma (mm)')
xlabel('Distance (cm)');
title('Spot size');
grid on

subplot(3,1,2);
errorbar(Zpos, measDose_Gy, dmeasDose_Gy,'r.');  hold on
errorbar(Zpos, predictedDose_Gy, dpredictedDose_Gy,'b.');  hold on
legend('Measured', 'Predicted');
ylabel('Mean dose (Gy)')
xlabel('Distance (cm)');
grid on
title('Central dose');

subplot(3,1,3);
errorbar(Zpos, ratios, dratios,'k.');  hold on
ylabel('Ratio')
xlabel('Distance (cm)');
grid on
title('Ratio measured / predicted');
ylim([0 1.2]);

%% Plottear resultados en función de LET
figure(3)
hold on
errorbar(Sw_preFilm, ratios, dratios, dratios, dSw_preFilm, dSw_preFilm, 'k.');
grid on
xlabel('LET (MeV/cm)')
ylabel('RE');
ylim([0 1.2]);
end

