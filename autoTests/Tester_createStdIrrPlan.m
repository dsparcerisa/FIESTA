clear all; close all

%% Plate
plateDose = nan(8, 12);
plateDose(2,:) = [nan nan nan nan 2 nan nan 8  nan nan nan nan];
plateDose(5,:) = [nan   1 nan nan 4 nan nan 12 nan nan nan nan];
NX = 12; NY = 8;
well2wellDist_cm = 0.899;
showPlate(plateDose)
title('CONV');

%% Microcubeta
% plateDose = [nan 1; 2 4; 6 8; 10 12];
% showPlate(plateDose)
% NX = 2; NY = 4;
% well2wellDist_cm = 1.125;

%% Create dose slice
E0 = 3;
z = 10;
dz = 0.001;
dxy = 0.01;
Nprot = 6.25e6; % 1 pC
targetTh = 0.001; % 10 um = 10^-5 m = 10^-3 cm
targetSPR = 1;
sizeX = 30;
sizeY = 30;
sigmaX = 0.1;
sigmaY = 0.1;
% Partimos de un haz con sigma = 0.1 mm
N0 = createGaussProfile(dxy, dxy, sizeX, sizeY, sigmaX, sigmaY);
doseSlice = getDoseMap(E0, z, dz, Nprot, targetTh, targetSPR, N0);

%% Find real sigma of doseSlice
Xsum = sum(doseSlice.data, 2);
F = fit((doseSlice.getAxisValues('X'))', Xsum, 'gauss1', 'StartPoint', [max(Xsum) 0 1]);
doseSigma = F.c1 / sqrt(2);

%% Calculate dose
deltaXY = doseSigma * 2.1;
I_FC1 = 1; % nA
I_factor = 0.81;
PP_factor = (10/250)^2; % PRUEBA
I_muestra = I_FC1 * I_factor * PP_factor; % nA

[plan, totIrrTime, doseRate] = createStdIrrPlan( plateDose, doseSlice, I_muestra, deltaXY, 4 );
plan.mode = 'CONV';
plan.Z = z*ones(size(plan.X))
fprintf('Total irradiation time in min: %f\n', totIrrTime/60);
fprintf('Dose rate: %f Gy/s\n', doseRate);

%%
plan.E = E0;
plan.Z = z;
scatter(plan.X(:), plan.Y(:), 100, plan.Q(:))
set(gca, 'XDir', 'reverse', 'YDir', 'reverse');

dose = getDoseFromPlan(CartesianGrid2D(N0), plan, dz, targetTh, targetSPR, N0);
figure;
dose.plotSlice
set(gca, 'Ydir', 'reverse', 'Xdir', 'reverse')

%% Calculate mean doses in all wells 
wellDiam = 6.35; % mm
wellRadius_cm = 0.1 * wellDiam / 2;

Nwells = sum(~isnan(plateDose(:)));
wells = {};

% Positions in reference with the center of the first spot
Xpos = well2wellDist_cm*(0:(-1):(-(NX-1)));
Ypos = well2wellDist_cm*(0:(NY-1));
[x,y] = meshgrid(Xpos, Ypos);

Xwells = x(~isnan(plateDose(:)));
Ywells = y(~isnan(plateDose(:)));

meanWellDoses = nan(Nwells, 1);
stdWellDoses = nan(Nwells, 1);

for i=1:Nwells
    well = getWell(CartesianGrid2D(dose), wellRadius_cm, [Xwells(i) Ywells(i)]);
    wellDoses = getStats(well, dose);
    meanWellDoses(i) = mean(wellDoses);   
    stdWellDoses(i) = std(wellDoses);   
end

meanWellDoses
stdWellDoses

%% Estimate what the EBT3 will look like
dose.data(isnan(dose.data))=0;
simRC = simulateRC(flip((dose.data)'+randn(size((dose.data)'))*0.3,2));
imshow(simRC)