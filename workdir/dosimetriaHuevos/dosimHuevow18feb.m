Z = 11;
inEgg_cm = 2.5;
Z_tumor =  Z + inEgg_cm;

%% Estimar la sigma sobre el tumor sin difusor

Z_18feb = 5:1.5:12.5;
sigma_18feb_noDif_mm = [0.944
1.513
1.851
2.437
2.856
3.011];
F18feb = fit(Z_18feb', sigma_18feb_noDif_mm, 'poly2')
sigma_noDif_tumor_mm = F18feb(Z_tumor)

%% Estimar el eefecto del difusor:

Z_19feb = 5.1:1.5:12.6;
sigma_18feb_noDif_mm = [0.478
0.644
0.915
1.148
1.267
1.431];
sigma_18feb_Dif_mm = [0.919
1.317
1.663
2.066
2.759
2.773];
effDif = sqrt(sigma_18feb_Dif_mm.^2 - sigma_18feb_noDif_mm.^2);

Fdif = fit(Z_19feb', effDif, 'poly2')
sigma_Dif_tumor_mm = Fdif(Z_tumor)

sigma_total_mm = sqrt(sigma_Dif_tumor_mm.^2 + sigma_noDif_tumor_mm.^2)

sigma_total_cm = sigma_total_mm / 10;

%% Creamos el cuadrado
pitch = 0.48; % cm

dxy = 0.01; % Tenth of mm
sizeXY = 8;

beamProfile = createGaussProfile(dxy, dxy, sizeXY, sizeXY, sigma_total_cm, sigma_total_cm);
doseCanvas = beamProfile.copy;
doseCanvas.data(:) = 0.0;
limits = [doseCanvas.minX doseCanvas.maxX doseCanvas.minY doseCanvas.maxY];

[X,Y] = meshgrid([-pitch/2 pitch/2]);
hold off
for i=1:4
    Di = beamProfile.copy;
    Di.shift([X(i) Y(i)]);
    doseCanvas = doseCanvas + Di;
    %doseCanvas.plotSlice
    %caxis([0 3e-4])
    %axis([-2 2 -2 2]);
    %colorbar
%     subplot(1,2,1)
%     hold on
%     plot(Di.getAxisValues('Y'),sum(Di.data),'b')
%     hold on
%     plot(doseCanvas.getAxisValues('Y'),sum(doseCanvas.data),'r');
%     subplot(1,2,2)
%     hold on
%     plot(Di.getAxisValues('X'),sum(Di.data,2),'b');
%     hold on
%     plot(doseCanvas.getAxisValues('X'),sum(doseCanvas.data,2),'r');    
%     pause
end
figure
doseCanvas.crop(limits);
doseCanvas.plotSlice

%% Los valores de beamProfile están normalizados de forma que la suma sea 1.

intensidadNominal = 0.160; % nA
factorConversion = 1; %medido con la FC
protonsXnA = 1e-9 / 1.602176634e-19;
timePerShot = 30; % seconds

totalProtons = doseCanvas.copy;
totalProtons.data = totalProtons.data * (intensidadNominal * factorConversion * protonsXnA * timePerShot);
totalProtons.plotSlice
axis([-0.25 0.25 -0.25 0.25]);

flux = totalProtons.copy;
flux.data = totalProtons.data / (dxy^2);
flux.plotSlice
axis([-0.5 0.5 -0.5 0.5]);
colorbar

%% Radius of tumours?
tumourMass = 0.2; % g;
tumourDensity = 1; %g/cm3 
tumourVolume = tumourMass/tumourDensity % cm3

% supposing spherical
tumourR_sph = (3 / 4 / pi * tumourVolume)^(1/3)
wellMap = getWell(totalProtons, tumourR_sph, [0 0]);
hold on;
contour(wellMap.getAxisValues('X'), wellMap.getAxisValues('Y'), wellMap.data)
title('Proton fluence (cm^-^2)');

figure(2)
allFluxes = flux.data(wellMap.data==1);
histogram(allFluxes)
meanVal = mean(allFluxes)
stdVal = std(allFluxes)
maxVal = max(allFluxes)

%% Sabiendo que el difusor de cobre mide 20 um de grosor, está a 27 mm de la salida, y se parte de una energía de 7.94 MeV a la salida de la ventana:
energyA = energyStoppingPowerCopper(7.94, 0:0.001:Z_tumor, 2.7, 20);
EA = energyA(end);

%% Energy loss por proton en los primeros 100 um
firstLayerThickness_cm = 0.01;
energyW = energyStoppingPowerWater(EA, 0:0.000001:firstLayerThickness_cm);
depE_MeV = EA - energyW(end)
depE_J = depE_MeV * 1e6 * 1.602176634e-19;

%% Por tanto el total de dosis es de
% D = flux * dE / rho * th
% Asumiendo densidad de agua
dose = flux.copy;
dose.data = flux.data * depE_J / (firstLayerThickness_cm * 0.001);
figure(2)
dose.plotSlice; colorbar
axis([-0.5 0.5 -0.5 0.5]);
max(dose.data(:))
figure(3)
allDoses = dose.data(wellMap.data==1);
histogram(allDoses)
meanVal = mean(allDoses)
stdVal = std(allDoses)
maxVal = max(allDoses)

%% Dosis media en el tumor
totProtInTumour = sum(totalProtons.data(wellMap.data==1))
totalEInTumour_J = totProtInTumour * EA * 1e6 * 1.602176634e-19;
totalTumourMass_kg = 4/3 * pi * tumourR_sph.^3 * 1e-3
meanTumourDose = totalEInTumour_J ./ totalTumourMass_kg