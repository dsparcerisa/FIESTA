polyX = [-6.132e-05 0.1326 -0.2029];
polyY = [-3.218e-06 0.1263 -0.07026];
Z = 7.5;
sigmaX = polyval(polyX, Z) / 10;
sigmaY = polyval(polyY, Z) / 10;
pitch = 0.12; % cm, 1.2 mm

dxy = 0.01; % Tenth of mm
sizeXY = 4;

beamProfile = createGaussProfile(dxy, dxy, sizeXY, sizeXY, sigmaX, sigmaY);
doseCanvas = beamProfile.copy;
doseCanvas.data(:) = 0.0;
limits = [doseCanvas.minX doseCanvas.maxX doseCanvas.minY doseCanvas.maxY];

[X,Y] = meshgrid([-0.12 0 0.12]);

for i=1:9

    Di = beamProfile.copy;
    Di.shift([X(i) Y(i)]);
    doseCanvas = doseCanvas + Di;
end

doseCanvas.crop(limits);
doseCanvas.plotSlice
axis([-2 2 -2 2]);
    
% Los valores de beamProfile están normalizados de forma que la suma sea 1.

intensidadNominal = 6; % nA
factorConversion = 1; %medido con la FC
protonsXnA = 1e-9 / 1.602176634e-19;
timePerShot = 10; % seconds

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
% supposing semipherical
tumourR_semisph = (1.5 / pi * tumourVolume)^(1/3)

% %% --> Asumimos una zona de 3x3 mm
% flux.crop([-0.3 0.3 -0.3 0.3]);
% allFluxes = flux.data(:);
% histogram(allFluxes)

%% Mejor asumimos un radio de 3mm
wellMap = getWell(totalProtons, tumourR_sph, [0 0]);
hold on;
contour(wellMap.getAxisValues('X'), wellMap.getAxisValues('Y'), wellMap.data, 'y')
title('Proton fluence (cm^-^2)');
%%
figure(2)
allFluxes = flux.data(wellMap.data==1);
histogram(allFluxes)
meanVal = mean(allFluxes)
stdVal = std(allFluxes)
maxVal = max(allFluxes)

%% Pérdida de energía en el kapton
E0 = 8; % MeV
Srho_kapton = 49.87; % MeVcm2 / g
rho_kapton = 1.43;
S_kapton = Srho_kapton * rho_kapton;
z_kapton = 8e-4; % cm
dE_kapton = z_kapton * S_kapton
E0_2 = E0 - dE_kapton

%% And now we consider deposited energy
[EA,~,stpW] = energyStoppingPower(E0_2,0:0.1:Z);
EA = EA(end)
stpW = stpW(end)

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
dose.plotSlice; colorbar
max(dose.data(:))

%% Pero la dosis media en el tumor es de :
totProtInTumour = sum(totalProtons.data(wellMap.data==1))
totalEInTumour_J = totProtInTumour * EA * 1e6 * 1.602176634e-19;
totalTumourMass_kg = 4/3 * pi * tumourR_sph.^3 * 1e-3;
meanTumourDose = totalEInTumour_J / totalTumourMass_kg

%% Valores de dosis esperados:
sigmaX = [
0.478
0.644
0.915
1.148
1.267
1.431];
sigmaY = [0.587
0.741
0.977
1.209
1.364
1.519];
sigmaX_cm = sigmaX / 10;
sigmaY_cm = sigmaY / 10;
Zval = [
5.1
6.6
8.1
9.6
11.1
12.6];
expDose = nan(6,1);
sigmaDose = nan(6,1);
EAval = nan(6,1);
Edepval = nan(6,1);
LET = nan(6,1);

totalNProt = 9.7 * 0.004 * 0.75 * protonsXnA;

% Capas de la RC
poliesterThk_cm = 0.0125; % cover
activeLayerThk_cm = 0.0028; % capa activa

innerRadius = 0.03; % 0.3 mm = 0.03 cm
wellMap = getWell(flux, innerRadius, [0 0]);

for i=1:6
    totalFluence_i = createGaussProfile(dxy, dxy, sizeXY, sizeXY, sigmaX_cm(i), sigmaY_cm(i));    
    totalFluence_i.data = totalNProt / dxy / dxy * totalFluence_i.data;
    
    % Degradar en aire
    [EA,~,stpW] = energyStoppingPower(E0_2,0:0.01:Zval(i));
    EA = EA(end)
    EAval(i) = EA;
    
    % Degradar en "agua", poliéster
    energyAL = energyStoppingPowerWater(EA, 0:0.000001:poliesterThk_cm);

    % Degradar en capa activa
    [energyW, Sw] = energyStoppingPowerWater(energyAL(end), 0:0.000001:activeLayerThk_cm);        
    depE_MeV = energyAL(end) - energyW(end)
    depE_J = depE_MeV * 1e6 * 1.602176634e-19;
    Edepval(i) = depE_MeV;
    LET(i) = Sw(end);
    
    % Calcular la dosis
    dose = totalFluence_i.copy;
    dose.data = totalFluence_i.data * depE_J / (activeLayerThk_cm * 0.001);
    
    % Calcular la dosis promedio en un radio dado
    allDoses = dose.data(wellMap.data==1);
    %histogram(allDoses);
    meanVal = mean(allDoses);
    stdVal = std(allDoses);
    maxVal = max(allDoses);
    expDose(i) = meanVal;
    sigmaDose(i) = stdVal;
end

%% Compared with measured data
measDose = [
92.504
60.459
31.543
20.804
14.934
11.677];
delta_measDose = [
8.338
15.700
2.243
1.417
1.052
0.834
];

errorbar(Zval, expDose, sigmaDose, 'r-')
hold on
errorbar(Zval, measDose, delta_measDose, 'b-')
NF =  median(measDose./expDose);
%errorbar(Zval, NF*expDose, NF*sigmaDose, 'm-')
legend('Calculated', 'Measured')%, 'Calculated-corrected');
grid on
ylabel('Dose (Gy)');
xlabel('Air depth (cm)');

%% Calcular eficiencia relativa
% Ref: Grilj, 2018
% Dhalf = a * LET + Dhalf_0
% RE = Dhalf_0 / Dhalf
% RE = Dhalf_o / (a * LET + Dhalf_0)
a = 1.31; % en Gy / (keV/um)
Dhalf_0 = 103; % Gy
RE = @(LET) Dhalf_0 ./ (a .* LET + Dhalf_0);
myRE = RE(LET/10);
relRE = myRE / myRE(1)

errorbar(Zval, expDose, sigmaDose, 'r-')
hold on
errorbar(Zval, measDose, delta_measDose, 'b-')
NF =  1./relRE;
errorbar(Zval, NF.*measDose, NF.*delta_measDose, 'm-')
legend('Calculated', 'Measured', 'Measured-corrected');
grid on
ylabel('Dose (Gy)');
xlabel('Air depth (cm)');