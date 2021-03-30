% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments

%% Abrir la GUI
stageControlStart(COMStage);

%% Alinear en posición 0 y medir la distancia (ANOTAR!)

beamExit2RCDistance = input('At Z=0, measure distance between beam exit and RC (cm): ');

%% Create dose plate

% Posible mejora: para tiempos del orden de <50ms, el tiempo de caida puede ser un issue. Mirar. 
plateDose = nan(8, 12);
plateDose(2,:) = [nan 1 nan nan 2 nan nan 3 nan nan 4 nan];
plateDose(5,:) = [nan 5 nan nan 6 nan nan 7 nan nan 8 nan];
plateDose(8,:) = [nan 9 nan nan 10 nan nan 11 nan nan 12 nan];
showPlate(plateDose)

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

%% Create plan
deltaXY = doseSigma * 2;
I_FC1 = 0.050; % 50 pA
I_factor = 0.76;
PP_factor = 1; % (10/250)^2; % PRUEBA
I_muestra = I_FC1 * I_factor * PP_factor; % nA

[plan, totIrrTime, doseRate] = createStdIrrPlan( plateDose, doseSlice, I_muestra, deltaXY, 9 );
fprintf('Total irradiation time in min: %f\n', totIrrTime/60);
fprintf('Dose rate: %f Gy/s\n', doseRate);

%% Calculate dose
plan.name = 'Plan de calibración';
plan.mode = 'CONV';
plan.codFiltro = '1';
plan.E = E0;
plan.Z = -(z-beamExit2RCDistance)*ones(size(plan.X));
plan.I = I_FC1;
scatter(plan.X(:), plan.Y(:), 100, plan.Q(:))
set(gca, 'XDir', 'reverse', 'YDir', 'reverse');

dose = getDoseFromPlan(CartesianGrid2D(N0), plan, dz, targetTh, targetSPR, N0, I_factor, beamExit2RCDistance);
figure;
dose.plotSlice
set(gca, 'Ydir', 'reverse', 'Xdir', 'reverse')

%% Save plan
writePlan(plan, 'plans/plan_calibrationRC.txt');

%% Calculate mean doses in all wells 
well2wellDist_cm = 0.899;
wellDiam = 6.35; % mm
wellRadius_cm = 0.1 * wellDiam / 2;

Nwells = sum(~isnan(plateDose(:)));
wells = {};

% Positions in reference with the center of the first spot
Xpos = well2wellDist_cm*(0:(-1):(-11));
Ypos = well2wellDist_cm*(0:(-1):(-7));
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
simRC = simulateRC(dose.data'+randn(size(dose.data'))*0.3);
imshow(simRC)